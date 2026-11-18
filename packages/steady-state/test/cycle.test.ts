import { describe, expect, it } from 'vitest'

import { F, Polynomial, solveCycle } from '../src/index.js'
import { bondGraph, enzyme, glut2, solveLinear } from './models.js'

const v = (name: string) => Polynomial.variable(name)

/** Random positive values for every variable a polynomial set uses. */
function randomValues(...ps: Polynomial[]): Record<string, number> {
  const values: Record<string, number> = {}
  for (const p of ps) for (const name of p.variables()) values[name] ??= 0.2 + 2 * Math.random()
  return values
}

describe('solveCycle', () => {
  it('derives the reversible Michaelis–Menten flux (Hunter et al. 2025, eq 8)', () => {
    const result = solveCycle(enzyme(), { cycleStates: ['q_E', 'q_C'], flux: 'v_R1', total: 'E0' })
    expect(result.ok).toBe(true)
    if (!result.ok) return
    const { flux } = result.value
    // κ1, κ2 the rates; K1 = K_S, K2 = K_P (solutes), K3 = K_E, K4 = K_C (enzyme)
    const [k1, k2, K1, K2, K3, K4, q1, q2] = ['kappa_R1', 'kappa_R2', 'K_S', 'K_P', 'K_E', 'K_C', 'q_S', 'q_P'].map(v) as Polynomial[]
    const num = v('E0').mul(k1!).mul(k2!).mul(K3!).mul(K4!).mul(K1!.mul(q1!).sub(K2!.mul(q2!)))
    const den = K4!.mul(k1!.add(k2!)).add(k1!.mul(K1!).mul(K3!).mul(q1!)).add(k2!.mul(K2!).mul(K3!).mul(q2!))
    expect(flux.num.mul(den).equals(num.mul(flux.den))).toBe(true)
    // The states' weights are the King–Altman ones: E ∝ κ1K4 + κ2K4, C ∝ κ1K1K3q1 + κ2K2K3q2
    const weights = result.value.weights
    const ratio = F.fraction(weights.get('q_C')!, weights.get('q_E')!)
    const expected = F.fraction(k1!.mul(K1!).mul(K3!).mul(q1!).add(k2!.mul(K2!).mul(K3!).mul(q2!)), K4!.mul(k1!.add(k2!)))
    expect(ratio.num.mul(expected.den).equals(expected.num.mul(ratio.den))).toBe(true)
  })

  it('gives the same flux whichever reaction of the cycle is reported', () => {
    const r1 = solveCycle(enzyme(), { cycleStates: ['q_E', 'q_C'], flux: 'v_R1' })
    const r2 = solveCycle(enzyme(), { cycleStates: ['q_E', 'q_C'], flux: 'v_R2' })
    expect(r1.ok && r2.ok).toBe(true)
    if (!r1.ok || !r2.ok) return
    expect(r1.value.flux.num.mul(r2.value.flux.den).equals(r2.value.flux.num.mul(r1.value.flux.den))).toBe(true)
  })

  it('solves the full GLUT2 cycle, agreeing with a numeric solve', () => {
    const states = ['q_1', 'q_2', 'q_3', 'q_4']
    const result = solveCycle(glut2(), { cycleStates: states, flux: 'v_R3' })
    expect(result.ok).toBe(true)
    if (!result.ok) return
    const { flux, matrix } = result.value
    for (let trial = 0; trial < 5; trial++) {
      const values = randomValues(flux.num, flux.den, ...matrix.flat())
      // Numerically: A q = 0 with the last row replaced by Σq = q_tot
      const A = matrix.map((row) => row.map((a) => a.evaluate(values)))
      A[3] = [1, 1, 1, 1]
      const q = solveLinear(A, [0, 0, 0, values.q_tot!])
      const numeric = values.kappa_R3! * (values.K_Ao! * values.q_Ao! * values.K_1! * q[0]! - values.K_2! * q[1]!)
      expect(F.evaluate(flux, values)).toBeCloseTo(numeric, 10)
    }
  })

  it('handles stoichiometry and slippage branches', () => {
    // A 2Na⁺ cotransporter-like cycle with a slippage reaction R5 between states 2 and 3
    const equations = bondGraph([
      { name: 'R1', forward: { E1: 1, Nao: 2 }, reverse: { E2: 1 } },
      { name: 'R2', forward: { E2: 1 }, reverse: { E3: 1 } },
      { name: 'R3', forward: { E3: 1 }, reverse: { E4: 1, Nai: 2 } },
      { name: 'R4', forward: { E4: 1 }, reverse: { E1: 1 } },
      { name: 'R5', forward: { E1: 1 }, reverse: { E4: 1 } },
    ])
    const states = ['q_E1', 'q_E2', 'q_E3', 'q_E4']
    const result = solveCycle(equations, { cycleStates: states, flux: 'v_R1' })
    expect(result.ok).toBe(true)
    if (!result.ok) return
    const { flux, matrix } = result.value
    expect(flux.num.degree('q_Nao')).toBe(2)
    const values = randomValues(flux.num, flux.den, ...matrix.flat())
    const A = matrix.map((row) => row.map((a) => a.evaluate(values)))
    A[3] = [1, 1, 1, 1]
    const q = solveLinear(A, [0, 0, 0, values.q_tot!])
    const numeric = values.kappa_R1! * (values.K_E1! * q[0]! * (values.K_Nao! * values.q_Nao!) ** 2 - values.K_E2! * q[1]!)
    expect(F.evaluate(flux, values)).toBeCloseTo(numeric, 10)
  })

  it('reports models it can’t reduce', () => {
    const code = (r: ReturnType<typeof solveCycle>) => (r.ok ? undefined : r.issues[0]!.code)
    // Dimerisation, 2E ⇌ D: the plain sum of the states isn't conserved
    const dimer = bondGraph([{ name: 'R1', forward: { E: 2 }, reverse: { D: 1 } }])
    expect(code(solveCycle(dimer, { cycleStates: ['q_E', 'q_D'], flux: 'v_R1' }))).toBe('no-conservation')
    const pair = bondGraph([
      { name: 'R1', forward: { E: 1, F: 1 }, reverse: { C: 1 } },
      { name: 'R2', forward: { C: 1 }, reverse: { E: 1, F: 1 } },
    ])
    expect(code(solveCycle(pair, { cycleStates: ['q_E', 'q_F', 'q_C'], flux: 'v_R1' }))).toBe('no-conservation')
    expect(code(solveCycle(enzyme(), { cycleStates: ['q_E', 'q_S'], flux: 'v_R1' }))).toBe('no-conservation')
    expect(code(solveCycle(enzyme(), { cycleStates: ['q_E', 'q_X'], flux: 'v_R1' }))).toBe('unknown-variable')
    // Two species that bind each other: no plain-sum conservation either
    const binding = bondGraph([
      { name: 'R1', forward: { A: 1, B: 1 }, reverse: { AB: 1 } },
      { name: 'R2', forward: { AB: 1 }, reverse: { A: 1, B: 1 } },
    ])
    expect(code(solveCycle(binding, { cycleStates: ['q_A', 'q_B', 'q_AB'], flux: 'v_R1' }))).toBe('no-conservation')
  })

  it('reports a term that isn’t linear in the states', () => {
    // E ⇌ F, E + F ⇌ G + E (catalysed by E itself), G ⇌ E: the total is conserved but R2 is quadratic
    const equations = bondGraph([
      { name: 'R1', forward: { E: 1 }, reverse: { F: 1 } },
      { name: 'R2', forward: { E: 1, F: 1 }, reverse: { G: 1, E: 1 } },
      { name: 'R3', forward: { G: 1 }, reverse: { E: 1 } },
    ])
    const r = solveCycle(equations, { cycleStates: ['q_E', 'q_F', 'q_G'], flux: 'v_R1' })
    expect(r.ok).toBe(false)
    if (!r.ok) expect(r.issues[0]!.code).toBe('nonlinear-in-states')
  })
})
