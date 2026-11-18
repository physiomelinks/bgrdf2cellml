import { cn, minus, plus, ci } from '@physiomelinks/math-ast'
import { describe, expect, it } from 'vitest'

import { Q, deriveSteadyState, evaluateAst, verify, type SteadyStateSpec } from '../src/index.js'
import { glut2, kir, nka, sglt1 } from './models.js'

const RT = 2578.73
const F = 96485

/** Evaluate an expression with the common form's variables and lumped parameters defined. */
function evaluateWith(result: ReturnType<typeof deriveSteadyState>, node: Parameters<typeof evaluateAst>[0], values: Record<string, number>) {
  const c = result.common!
  return evaluateAst(node, values, new Map([...c.variables, ...c.definitions].map(({ name, value }) => [name, value])))
}

describe('the common form', () => {
  const glut2Spec: SteadyStateSpec = { label: 'GLUT2', cycleStates: ['q_1', 'q_2', 'q_3', 'q_4'], flux: 'v_R2', fast: ['kappa_R3', 'kappa_R4'] }

  it('writes GLUT2 as κ/Φ (q̄o − q̄i) with Φ = 1 + q̄o/k + q̄i/k + q̄o q̄i/k (eqs 19–20)', () => {
    const result = deriveSteadyState(glut2(), glut2Spec)
    expect(result.issues).toEqual([])
    const c = result.common!
    expect(c.bracket).toBe('difference')
    expect(c.qeq).toBeUndefined()
    expect(c.exponent).toBeUndefined()
    expect(c.reactants).toEqual([['qbar_Ao', 1]])
    expect(c.products).toEqual([['qbar_Ai', 1]])
    expect(c.phi.terms.map((t) => t.k?.name)).toEqual(expect.arrayContaining(['k_Ao_GLUT2', 'k_Ai_GLUT2', 'k_Ai_Ao_GLUT2']))
    // Table 3's parameters give eq 20's km values
    const values = { K_1: 33.2, K_2: 4.25e3, K_3: 344.59, K_4: 1.99, kappa_R1: 0.36, kappa_R2: 0.26, q_tot: 1 }
    const k = (name: string) => evaluateWith(result, ci(name), values)
    expect(k('k_Ao_GLUT2')).toBeCloseTo(1.4667, 3)
    expect(k('k_Ai_GLUT2')).toBeCloseTo(21.603, 2)
    expect(k('k_Ai_Ao_GLUT2')).toBeCloseTo(235.08, 1)
    expect(verify(glut2(), glut2Spec, result, { fastScale: 1e8 }).maxRelativeError).toBeLessThan(1e-5)
    expect(verify(glut2(), glut2Spec, result).commonFormError).toBeLessThan(1e-12)
  })

  it('writes GLUT2 in concentrations (sodium homeostasis paper, eq 4 and Table 1)', () => {
    const spec = { ...glut2Spec, concentrations: { q_Ao: 'V_o', q_Ai: 'V_i' } }
    const result = deriveSteadyState(glut2(), spec)
    const values = { K_1: 33.2, K_2: 4.25e3, K_3: 344.59, K_4: 1.99, kappa_R1: 0.36, kappa_R2: 0.26, q_tot: 1, K_Ao: 149.65, K_Ai: 149.65, V_o: 0.09, V_i: 0.09 }
    expect(evaluateWith(result, ci('k_Ao_GLUT2'), values)).toBeCloseTo(0.1089, 3) // 1.4667/(149.65·0.09)
    expect(evaluateWith(result, ci('k_Ai_GLUT2'), values)).toBeCloseTo(1.604, 2)
    expect(result.common!.qeq).toBeDefined() // K_Ai V_i/(K_Ao V_o): 1 with these values
    expect(evaluateWith(result, ci('Q_eq_GLUT2'), values)).toBeCloseTo(1, 12)
    expect(verify(glut2(), spec, result).commonFormError).toBeLessThan(1e-12)
  })

  it('writes SGLT1 as a ratio with z = 2 and agrees with eqs 37–38', () => {
    const spec: SteadyStateSpec = { label: 'SGLT1', cycleStates: ['q_1', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6'], flux: 'v_R3', fast: ['kappa_R1', 'kappa_R2', 'kappa_R4', 'kappa_R5'] }
    const result = deriveSteadyState(sglt1(), spec)
    expect(result.issues).toEqual([])
    const c = result.common!
    expect(c.bracket).toBe('ratio')
    expect(c.exponent?.z).toEqual(Q.rational(2))
    expect(c.exponent?.potential).toBe('u_m')
    // Eq 37–38 with Table 5's parameters (z1 = z2 = 1/2)
    const p: Record<string, number> = { K_1: 2.235, K_2: 10.44, K_3: 8.602, K_4: 8.602, K_5: 47.71, K_6: 0.319, kappa_R3: 5.812, kappa_R6: 15.66, q_tot: 1, RT, F, u_m: -0.05 }
    const qb = { Nao: 2.1, Nai: 0.3, Glco: 0.8, Glci: 0.4 }
    const values = { ...p, K_Nao: 1, K_Nai: 1, K_Glco: 1, K_Glci: 1, q_Nao: qb.Nao, q_Nai: qb.Nai, q_Glco: qb.Glco, q_Glci: qb.Glci }
    const e = (x: number) => Math.exp((x * F * p.u_m!) / RT)
    const B =
      (e(1) + qb.Nao ** 2 * (p.K_1! / p.K_2! + (p.K_1! / p.K_3!) * qb.Glco)) * (qb.Nai ** 2 * qb.Glci * e(0.5) + p.kappa_R6! / p.kappa_R3!) +
      (p.K_1! / p.K_6! + qb.Nai ** 2 * (p.K_1! / p.K_5! + (p.K_1! / p.K_4!) * qb.Glci)) * (qb.Nao ** 2 * qb.Glco * e(0.5) + (p.kappa_R6! / p.kappa_R3!) * e(2))
    const eq37 = (p.q_tot! * p.kappa_R6! * p.K_1! * (qb.Nao ** 2 * qb.Glco - qb.Nai ** 2 * qb.Glci * e(2))) / B
    expect(evaluateWith(result, c.flux, values) / eq37).toBeCloseTo(1, 12)
    // The reversal potential, RT/(2F) ln((q̄Nao/q̄Nai)² q̄Glco/q̄Glci) (eq 9)
    const uRev = evaluateWith(result, c.reversalPotential!, values)
    expect(uRev).toBeCloseTo((RT / (2 * F)) * Math.log((qb.Nao / qb.Nai) ** 2 * (qb.Glco / qb.Glci)), 12)
    expect(verify(sglt1(), spec, result, { fastScale: 1e9 }).maxRelativeError).toBeLessThan(1e-5)
  })

  it('asks for held species when the model has none as states', () => {
    // Michaelis–Menten with S and P as constants
    const d = (q: string) => ({ type: 'Derivative' as const, expression: ci(q), variable: ci('t') })
    const v1 = { type: 'Multiply' as const, children: [ci('k1'), minus({ type: 'Multiply' as const, children: [ci('E'), ci('S')] }, ci('C'))] }
    const v2 = { type: 'Multiply' as const, children: [ci('k2'), minus(ci('C'), { type: 'Multiply' as const, children: [ci('E'), ci('P')] })] }
    const equations = [
      { type: 'Equal' as const, left: ci('v1'), right: v1 },
      { type: 'Equal' as const, left: ci('v2'), right: v2 },
      { type: 'Equal' as const, left: d('E'), right: minus(ci('v2'), ci('v1')) },
      { type: 'Equal' as const, left: d('C'), right: minus(ci('v1'), ci('v2')) },
    ]
    const without = deriveSteadyState(equations, { cycleStates: ['E', 'C'], flux: 'v1' })
    expect(without.issues.map((i) => i.code)).toEqual(['no-held-species'])
    const named = deriveSteadyState(equations, { cycleStates: ['E', 'C'], flux: 'v1', solutes: ['S', 'P'] })
    expect(named.issues).toEqual([])
    expect(named.common!.reactants).toEqual([['S', 1]])
  })

  it('reports SGLT1 with slippage as exact but not in the common form', () => {
    const spec: SteadyStateSpec = { cycleStates: ['q_1', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6'], flux: 'v_R3' }
    const result = deriveSteadyState(sglt1(0.5, 0.5, true), spec)
    expect(result.ok).toBe(true)
    expect(result.issues.map((i) => i.code)).toEqual(['not-thermodynamically-consistent'])
    expect(verify(sglt1(0.5, 0.5, true), spec, result).maxRelativeError).toBeLessThan(1e-9)
  })

  it('writes Kir with z = −1 once charge balance is given (eqs 19–23)', () => {
    const relations = { z4: plus(minus(cn(-1), ci('z2')), ci('z1'), ci('z3'), ci('z5'), minus(ci('z6'))) }
    const spec: SteadyStateSpec = { label: 'Kir', cycleStates: ['q_1', 'q_2', 'q_3'], flux: 'v_R2', fast: ['kappa_R1', 'kappa_R3'], opaque: ['u_m'], relations }
    const result = deriveSteadyState(kir(), spec)
    expect(result.issues).toEqual([])
    const c = result.common!
    expect(c.bracket).toBe('difference')
    expect(c.exponent?.z).toEqual(Q.rational(-1))
    const values = { z1: 0, z2: 0, z3: 0, z5: 0, z6: -0.82, K_1: 2, K_2: 0.7, K_3: 1.4, kappa_R2: 0.9, q_tot: 1, K_Ki: 0.0239, K_Ko: 0.0239, q_Ki: 140, q_Ko: 4.5, RT, F, u_m: -0.06 }
    // u_rev = −(RT/F) ln([K]i/[K]o) = −91.2 mV at RT/F = 26.7 mV
    expect(evaluateWith(result, c.reversalPotential!, values)).toBeCloseTo(-(RT / F) * Math.log(140 / 4.5), 12)
    expect(verify(kir(), spec, result, { fastScale: 1e9 }).maxRelativeError).toBeLessThan(1e-5)
    expect(verify(kir(), spec, result).commonFormError).toBeLessThan(1e-9)
  })

  it('writes NKA with ATP hydrolysis as the driver and gives u_rev = −225 mV (eqs 10–12)', () => {
    const spec: SteadyStateSpec = {
      label: 'NKA',
      cycleStates: ['q_1', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6'],
      flux: 'v_R6',
      driver: ['q_ATP', 'q_H2O', 'q_ADP', 'q_Pi', 'q_H'],
      relations: { z2: minus(cn(1), ci('z1')) },
    }
    const result = deriveSteadyState(nka(), spec)
    expect(result.issues).toEqual([])
    const c = result.common!
    expect(c.bracket).toBe('ratio')
    expect(c.exponent?.z).toEqual(Q.rational(-1))
    expect(c.driver?.reactants.map(([v]) => v)).toEqual(['qbar_ATP', 'qbar_H2O'])
    // Q_ATP/Q_ATP^eq = 2.3e10, [Na]i/[Na]o = 10/140, [K]o/[K]i = 4.5/145, RT/F = 25 mV
    const RTp = 25e-3 * F
    const values: Record<string, number> = { RT: RTp, F, z1: 0.945, u_m: -0.06, q_tot: 1 }
    for (const s of ['Nai', 'Nao', 'Ki', 'Ko', 'ATP', 'H2O', 'ADP', 'Pi', 'H']) values[`K_${s}`] = 1
    Object.assign(values, { q_Nai: 10, q_Nao: 140, q_Ki: 145, q_Ko: 4.5, q_ATP: 2.3e10, q_H2O: 1, q_ADP: 1, q_Pi: 1, q_H: 1 })
    for (const [i, K] of [0.01673, 4.477e3, 46.64e3, 720.6, 1.607, 30.25e3].entries()) values[`K_${i + 1}`] = K
    for (const [i, k] of [2.839e3, 51.1, 469.3, 619.1e3, 35.87, 868.4].entries()) values[`kappa_R${i + 1}`] = k
    expect(evaluateWith(result, c.reversalPotential!, values) * 1e3).toBeCloseTo(-224.9, 1)
    expect(verify(nka(), spec, result).maxRelativeError).toBeLessThan(1e-9)
    expect(verify(nka(), spec, result).commonFormError).toBeLessThan(1e-9)
  })
})
