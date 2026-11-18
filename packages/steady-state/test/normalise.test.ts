import { ci, cn, divide, exp, fn, minus, plus, power, times, type AstNode } from '@physiomelinks/math-ast'
import { describe, expect, it } from 'vitest'

import { F, Polynomial, SteadyStateError, atomValues, evaluateAst, fractionToAst, normalise } from '../src/index.js'

const ln = (a: AstNode) => fn('ln', a)
/** `u = RT ln(K q)`, as bond-graph CellML writes a chemical potential */
const potential = (K: string, q: string): AstNode => times(ci('RT'), ln(times(ci(K), ci(q))))
/** `κ (exp(u_f/RT) − exp(u_r/RT))` */
const rate = (kappa: string, uf: AstNode, ur: AstNode): AstNode =>
  times(ci(kappa), minus(exp(divide(uf, ci('RT'))), exp(divide(ur, ci('RT')))))

const v = (name: string) => Polynomial.variable(name)

/** Check that the normalised fraction and the original expression agree numerically. */
function agree(node: AstNode, definitions: Map<string, AstNode>, values: Record<string, number>, states: string[] = [], opaque: string[] = []) {
  const { fractions, atoms } = normalise([node], { definitions, states, opaque })
  const f = fractions[0]!
  const expected = evaluateAst(node, values, definitions)
  const env = atomValues(atoms, { ...values, ...Object.fromEntries([...definitions.keys()].filter((k) => !(k in values)).map((k) => [k, evaluateAst(ci(k), values, definitions)])) })
  expect(F.evaluate(f, env)).toBeCloseTo(expected, 9)
  // and written back as an expression tree
  expect(evaluateAst(fractionToAst(f, atoms), values, definitions)).toBeCloseTo(expected, 9)
  return { f, atoms }
}

describe('normalise', () => {
  it('turns a bond-graph reaction into mass action, cancelling RT', () => {
    const definitions = new Map<string, AstNode>([
      ['u_E', potential('K_E', 'q_E')],
      ['u_S', potential('K_S', 'q_S')],
      ['u_C', potential('K_C', 'q_C')],
      ['u_f', plus(ci('u_E'), ci('u_S'))],
      ['u_r', ci('u_C')],
      ['v_R1', rate('kappa', ci('u_f'), ci('u_r'))],
    ])
    const { fractions, atoms } = normalise([ci('v_R1')], { definitions, states: ['q_E', 'q_C'] })
    expect(atoms.size).toBe(0)
    const expected = v('kappa').mul(v('K_E').mul(v('q_E')).mul(v('K_S')).mul(v('q_S')).sub(v('K_C').mul(v('q_C'))))
    expect(F.isPolynomial(fractions[0]!)).toBe(true)
    expect(fractions[0]!.num.equals(expected)).toBe(true)
  })

  it('raises bond counts to powers', () => {
    const definitions = new Map<string, AstNode>([
      ['u_1', potential('K', 'q_1')],
      ['u_2', potential('K', 'q_2')],
    ])
    const node = exp(divide(plus(ci('u_1'), times(cn(2), ci('u_2'))), ci('RT')))
    const { fractions } = normalise([node], { definitions })
    expect(fractions[0]!.num.equals(v('K').pow(3).mul(v('q_1')).mul(v('q_2').pow(2)))).toBe(true)
    agree(node, definitions, { RT: 2578.73, K: 2, q_1: 0.3, q_2: 0.7 })
  })

  it('keeps a membrane potential in exp atoms and doesn’t inline it', () => {
    // The fork's electrogenic reaction: a transformer ratio on the membrane potential
    const definitions = new Map<string, AstNode>([
      ['u_C_m', divide(ci('q_C_m'), ci('C_C_m'))],
      ['u_q_B', potential('K_q_B', 'q_B')],
      ['u_TF_f', times(ci('TF_f'), ci('u_C_m'))],
      ['u_TF_r', divide(ci('u_C_m'), ci('TF_r'))],
      ['u_f', plus(ci('u_Se_A'), ci('u_TF_f'))],
      ['u_r', plus(ci('u_TF_r'), ci('u_q_B'))],
      ['v_Re', rate('kappa', ci('u_f'), ci('u_r'))],
    ])
    const { f, atoms } = agree(
      ci('v_Re'),
      definitions,
      { RT: 2578.73, kappa: 0.5, K_q_B: 3, q_B: 0.2, q_C_m: 1e-3, C_C_m: 1e-2, TF_f: -48242, TF_r: -2.07e-5, u_Se_A: 1200 },
      ['q_B'],
      ['u_C_m'],
    )
    // exp(u_Se_A/RT), exp(TF_f u_C_m/RT) and exp(u_C_m/(TF_r RT))
    expect([...atoms.values()].map((a) => a.base.map(([x, e]) => `${x}^${e}`).join(' '))).toEqual([
      'RT^-1 TF_f^1 u_C_m^1',
      'RT^-1 TF_r^-1 u_C_m^1',
      'RT^-1 u_Se_A^1',
    ])
    expect(f.num.variables()).not.toContain('q_C_m')
    // Not opaque, the exponents are written in the charge
    const { atoms: inCharge } = normalise([ci('v_Re')], { definitions, states: ['q_B'] })
    expect([...inCharge.values()].some((a) => a.base.some(([x]) => x === 'q_C_m'))).toBe(true)
  })

  it('shares one atom between multiples of the same exponent', () => {
    const half = exp(divide(times(cn(0.5), ci('F'), ci('u')), ci('RT')))
    const whole = exp(divide(times(ci('F'), ci('u')), ci('RT')))
    const { fractions, atoms } = normalise([times(half, whole), divide(power(half, cn(2)), whole)], { definitions: new Map() })
    expect(atoms.size).toBe(1)
    const [atom] = [...atoms.values()]
    expect(atom!.denominator).toBe(2)
    expect(fractions[0]!.num.equals(v(atom!.name).pow(3))).toBe(true)
    expect(fractions[1]!.num.equals(Polynomial.one)).toBe(true)
  })

  it('represents a symbolic charge as part of the atom', () => {
    // Kir: exp(z1 F u/RT)·exp(−z2 F u/RT), with z1, z2 unknown
    const node = times(exp(divide(times(ci('z1'), ci('F'), ci('u')), ci('RT'))), exp(minus(divide(times(ci('z2'), ci('F'), ci('u')), ci('RT')))))
    const { fractions, atoms } = normalise([node], { definitions: new Map() })
    expect(atoms.size).toBe(2)
    expect(F.isPolynomial(fractions[0]!)).toBe(false)
    agree(node, new Map(), { z1: 0.3, z2: 0.2, F: 96485, u: -0.06, RT: 2578.73 })
  })

  it('keeps fractions exact', () => {
    const node = divide(plus(ci('a'), ci('b')), times(cn(4), ci('a')))
    const { fractions } = normalise([plus(node, cn(0.1))], { definitions: new Map() })
    const f = fractions[0]!
    expect(f.den.equals(v('a').scale(20))).toBe(false) // the denominator's content is 1
    expect(f.den.equals(v('a'))).toBe(true)
    agree(plus(node, cn(0.1)), new Map(), { a: 0.7, b: 1.9 })
  })

  it('reports what it can’t represent', () => {
    const issue = (node: AstNode, definitions = new Map<string, AstNode>()) => {
      try {
        normalise([node], { definitions })
      } catch (e) {
        expect(e).toBeInstanceOf(SteadyStateError)
        return (e as SteadyStateError).issue.code
      }
      return undefined
    }
    expect(issue(fn('sin', ci('x')))).toBe('unsupported-construct')
    expect(issue(ln(ci('x')))).toBe('unsupported-logarithm')
    expect(issue(exp(times(ci('x'), ln(ci('y')))))).toBe('unsupported-logarithm')
    expect(issue(exp(times(cn(0.5), ln(ci('y')))))).toBe('unsupported-logarithm')
    expect(issue(power(ci('x'), ci('n')))).toBe('unsupported-construct')
    expect(issue(ci('a'), new Map<string, AstNode>([['a', plus(ci('b'), ln(ci('a')))]]))).toBe('circular-definition')
  })
})

describe('writing polynomials', () => {
  it('takes common factors out of a sum', async () => {
    const { polynomialToAst } = await import('../src/index.js')
    const [K1, K2, K4, k1] = ['K_1', 'K_2', 'K_4', 'kappa_1'].map(v)
    const p = K1!.mul(K2!).mul(k1!).add(K2!.mul(K4!).mul(k1!)).scale(2)
    const node = polynomialToAst(p)
    // 2 · K_2 · κ_1 · (K_1 + K_4)
    expect(node.type).toBe('Multiply')
    expect(node.type === 'Multiply' && node.children.at(-1)!.type).toBe('Add')
    expect(evaluateAst(node, { K_1: 2, K_2: 3, K_4: 5, kappa_1: 7 })).toBeCloseTo(2 * 3 * 7 * (2 + 5), 12)
    // A sum with nothing in common is left as it is
    expect(polynomialToAst(K1!.add(K4!)).type).toBe('Add')
  })
})
