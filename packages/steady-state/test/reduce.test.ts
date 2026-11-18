import { ci, cn, minus, plus, walk, type AstNode, type Equation } from '@physiomelinks/math-ast'
import { describe, expect, it } from 'vitest'

import { deriveSteadyState, evaluateAst, reduceModel, splitEquations, suggestSpecs, verify, type SteadyStateSpec } from '../src/index.js'
import { enzyme, glut2, kir, nka, sglt1 } from './models.js'

const name = (e: Equation) => (e.left.type === 'Identifier' ? e.left.name : e.left.type === 'Derivative' && e.left.expression.type === 'Identifier' ? `d${e.left.expression.name}` : '?')

describe('suggestSpecs', () => {
  it('finds GLUT2’s cycle, its held glucose and its binding steps', () => {
    const [cycle, ...rest] = suggestSpecs(glut2())
    expect(rest).toEqual([])
    expect([...cycle!.spec.cycleStates].sort()).toEqual(['q_1', 'q_2', 'q_3', 'q_4'])
    expect(cycle!.solutes.sort()).toEqual(['q_Ai', 'q_Ao'])
    expect(cycle!.spec.fast).toEqual(['kappa_R3', 'kappa_R4'])
    expect(cycle!.spec.flux).toBe('v_R1')
  })

  it('finds SGLT1’s six states and four binding steps, not a glucose moiety', () => {
    const [cycle] = suggestSpecs(sglt1())
    expect([...cycle!.spec.cycleStates].sort()).toEqual(['q_1', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6'])
    expect(cycle!.spec.fast).toEqual(['kappa_R1', 'kappa_R2', 'kappa_R4', 'kappa_R5'])
  })

  it('suggests no fast reactions when every step binds something (NKA, Michaelis–Menten)', () => {
    expect(suggestSpecs(nka())[0]!.spec.fast).toBeUndefined()
    const [mm] = suggestSpecs(enzyme())
    expect([...mm!.spec.cycleStates].sort()).toEqual(['q_C', 'q_E'])
    expect(mm!.spec.fast).toBeUndefined()
  })

  it('recognises transformer ratios as multiples of F', () => {
    const equations = glut2()
    const [cycle] = suggestSpecs(equations, { ratios: ['TF_f', 'TF_r'], values: { TF_f: -48242.666, TF_r: 2.07285049e-5, F: 96485.332 } })
    expect(cycle!.spec.relations).toBeDefined()
    expect(evaluateAst(cycle!.spec.relations!.TF_f!, { F: 96485.332 })).toBeCloseTo(-48242.666, 2)
    expect(evaluateAst(cycle!.spec.relations!.TF_r!, { F: 96485.332 })).toBeCloseTo(2.07285049e-5, 9)
  })
})

describe('reduceModel', () => {
  const values: Record<string, number> = { RT: 2578.73, q_tot: 1 }

  /** Each held species' rate in the reduced model equals the full model's at its steady state. */
  function checkHeld(full: Equation[], spec: SteadyStateSpec, extra: Record<string, number> = {}) {
    const result = deriveSteadyState(full, spec)
    expect(result.ok).toBe(true)
    const reduced = reduceModel(full, spec, result)
    const { definitions, rates } = splitEquations(reduced.equations as Equation[])
    // No cycle state's equation is left
    for (const q of spec.cycleStates) expect(rates.has(q)).toBe(false)
    // Evaluate with random values and compare with the steady-state flux
    const free = new Map<string, number>(Object.entries({ ...values, ...extra }))
    const valueOf = (node: AstNode) => {
      for (;;)
        try {
          return evaluateAst(node, free, definitions)
        } catch (e) {
          const missing = /No value for (\S+)/.exec(String(e))?.[1]
          if (missing === undefined) throw e
          free.set(missing, 0.5 + Math.random())
        }
    }
    const v = valueOf(ci(reduced.flux))
    expect(Number.isFinite(v)).toBe(true)
    return { reduced, rates, valueOf, v, result }
  }

  it('reduces GLUT2: glucose moves at ∓v and the transporter states are gone', () => {
    const spec: SteadyStateSpec = { label: 'GLUT2', cycleStates: ['q_1', 'q_2', 'q_3', 'q_4'], flux: 'v_R1', fast: ['kappa_R3', 'kappa_R4'] }
    const { reduced, rates, valueOf, v } = checkHeld(glut2(), spec)
    expect(reduced.flux).toBe('v_GLUT2')
    expect(valueOf(rates.get('q_Ao')!)).toBeCloseTo(-v, 12)
    expect(valueOf(rates.get('q_Ai')!)).toBeCloseTo(v, 12)
    expect(reduced.parameters).toEqual(['q_tot'])
    expect(reduced.sections.map((s) => s.title)).toEqual(['Steady-state flux', 'Held species'])
    // The cycle's own definitions (its potentials and rates) are gone; glucose's potentials only the cycle used
    const names = reduced.equations.map(name)
    for (const gone of ['u_1', 'v_R3', 'u_f_R1', 'u_Ao']) expect(names).not.toContain(gone)
  })

  it('reduces SGLT1: two sodium per glucose', () => {
    const spec: SteadyStateSpec = { label: 'SGLT1', cycleStates: ['q_1', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6'], flux: 'v_R3', fast: ['kappa_R1', 'kappa_R2', 'kappa_R4', 'kappa_R5'] }
    const { rates, valueOf, v } = checkHeld(sglt1(), spec, { F: 96485, u_m: -0.05 })
    expect(valueOf(rates.get('q_Nao')!)).toBeCloseTo(-2 * v, 12)
    expect(valueOf(rates.get('q_Nai')!)).toBeCloseTo(2 * v, 12)
    expect(valueOf(rates.get('q_Glci')!)).toBeCloseTo(v, 12)
  })

  it('reduces Kir with its charge relation', () => {
    const relations = { z4: plus(minus(cn(-1), ci('z2')), ci('z1'), ci('z3'), ci('z5'), minus(ci('z6'))) }
    const spec: SteadyStateSpec = { label: 'Kir', cycleStates: ['q_1', 'q_2', 'q_3'], flux: 'v_R2', fast: ['kappa_R1', 'kappa_R3'], opaque: ['u_m'], relations }
    const { rates, valueOf, v, result } = checkHeld(kir(), spec, { F: 96485, u_m: -0.06 })
    expect(valueOf(rates.get('q_Ko')!)).toBeCloseTo(v, 12)
    expect(verify(kir(), spec, result).commonFormError).toBeLessThan(1e-9)
  })

  it('keeps the rest of a model as it was', () => {
    // GLUT2 with a downstream reaction consuming inside glucose: Ai → B
    const extra: Equation[] = [
      { type: 'Equal', left: ci('v_use'), right: { type: 'Multiply', children: [ci('k_use'), ci('q_Ai')] } },
      { type: 'Equal', left: { type: 'Derivative', expression: ci('q_B'), variable: ci('t') }, right: ci('v_use') },
    ]
    const full = glut2().map((e) =>
      e.left.type === 'Derivative' && e.left.expression.type === 'Identifier' && e.left.expression.name === 'q_Ai' ? { ...e, right: minus(e.right, ci('v_use')) } : e,
    )
    const spec: SteadyStateSpec = { cycleStates: ['q_1', 'q_2', 'q_3', 'q_4'], flux: 'v_R1' }
    const { reduced, rates, valueOf, v } = checkHeld([...full, ...extra], spec)
    expect(reduced.equations.map(name)).toEqual(expect.arrayContaining(['v_use', 'dq_B']))
    expect(valueOf(rates.get('q_Ai')!)).toBeCloseTo(v - valueOf(ci('v_use')), 12)
  })
})

describe('reduceModel keep', () => {
  it('keeps a reaction’s flux at its steady state, so what names it still can', () => {
    const spec: SteadyStateSpec = { label: 'SGLT1', cycleStates: ['q_1', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6'], flux: 'v_R3', fast: ['kappa_R1', 'kappa_R2', 'kappa_R4', 'kappa_R5'] }
    const full = sglt1()
    const result = deriveSteadyState(full, spec)
    const reduced = reduceModel(full, spec, result, { keep: ['v_R1', 'v_R5'] })
    const { definitions, rates } = splitEquations(reduced.equations as Equation[])
    // Each kept flux is the flux itself (once per cycle), and the held species' rates use it
    expect(definitions.get('v_R1')).toEqual(ci('v_SGLT1'))
    expect(definitions.get('v_R5')).toEqual(ci('v_SGLT1'))
    expect([...walk(rates.get('q_Nao')!)].some((n) => n.type === 'Identifier' && n.name === 'v_R1')).toBe(true)
    expect(reduced.sections.map((s) => s.title)).toContain('Kept at steady state')
    // A kept variable isn't a parameter the model lacks
    expect(reduced.parameters).toEqual(['q_tot'])
  })
})
