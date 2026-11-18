import { ci, times } from '@physiomelinks/math-ast'
import { describe, expect, it } from 'vitest'

import { F, Polynomial, atomValues, solveCycle, type Fraction } from '../src/index.js'
import { bondGraph, glut2 } from './models.js'

const v = (name: string) => F.fromPolynomial(Polynomial.variable(name))
const c = (x: number) => F.fromPolynomial(Polynomial.constant(x))
const same = (a: Fraction, b: Fraction) => a.num.mul(b.den).equals(b.num.mul(a.den))

describe('rapid equilibrium', () => {
  it('reproduces GLUT2’s steady-state flux (Hunter et al. 2025, eqs 19–20)', () => {
    const states = ['q_1', 'q_2', 'q_3', 'q_4']
    for (const flux of ['v_R1', 'v_R2', 'v_R3', 'v_R4']) {
      const result = solveCycle(glut2(), { cycleStates: states, flux, fast: ['kappa_R3', 'kappa_R4'] })
      expect(result.ok).toBe(true)
      if (!result.ok) return
      const [K1, K2, K3, K4, k1, k2] = ['K_1', 'K_2', 'K_3', 'K_4', 'kappa_R1', 'kappa_R2'].map(v) as Fraction[]
      const qo = F.mul(v('K_Ao'), v('q_Ao'))
      const qi = F.mul(v('K_Ai'), v('q_Ai'))
      const L = F.div(k2!, k1!)
      const r = (a: Fraction, b: Fraction) => F.div(a, b)
      const den = [
        F.add(r(K3!, K1!), r(K3!, K4!)),
        F.mul(F.add(r(K3!, K2!), F.mul(L, r(K3!, K4!))), qo),
        F.mul(F.add(c(1), F.mul(L, r(K3!, K1!))), qi),
        F.mul(F.mul(L, F.add(c(1), r(K3!, K2!))), F.mul(qo, qi)),
      ].reduce(F.add)
      const expected = F.div(F.mul(F.mul(F.mul(v('q_tot'), k2!), K3!), F.sub(qo, qi)), den)
      expect(same(result.value.flux, expected)).toBe(true)
      // No fast rate constant survives
      expect(result.value.flux.num.has('kappa_R3') || result.value.flux.den.has('kappa_R4')).toBe(false)
    }
  })

  it('gives the paper’s fitted km values from Table 3', () => {
    const result = solveCycle(glut2(), { cycleStates: ['q_1', 'q_2', 'q_3', 'q_4'], flux: 'v_R2', fast: ['kappa_R3', 'kappa_R4'] })
    if (!result.ok) throw new Error(result.issues[0]!.message)
    const values = { K_1: 33.2, K_2: 4.25e3, K_3: 344.59, K_4: 1.99, kappa_R1: 0.36, kappa_R2: 0.26, q_tot: 1, K_Ao: 1, K_Ai: 1 }
    const flux = (qo: number, qi: number) => F.evaluate(result.value.flux, { ...values, q_Ao: qo, q_Ai: qi })
    // v(q̄o, 0) = vm·q̄o/(1 + q̄o/km¹): km¹ ≈ 1.4667 is where the flux is half its maximum
    const vmax = flux(1e9, 0)
    expect(flux(1.4667, 0) / vmax).toBeCloseTo(0.5, 3)
    expect(-flux(0, 21.603) / -flux(0, 1e9)).toBeCloseTo(0.5, 3)
  })

  it('agrees with the full cycle as the fast rate constants grow', () => {
    const states = ['q_1', 'q_2', 'q_3', 'q_4']
    const limit = solveCycle(glut2(), { cycleStates: states, flux: 'v_R2', fast: ['kappa_R3', 'kappa_R4'] })
    const full = solveCycle(glut2(), { cycleStates: states, flux: 'v_R2' })
    if (!limit.ok || !full.ok) throw new Error('not solved')
    const values: Record<string, number> = { K_1: 3, K_2: 0.5, K_3: 2, K_4: 1.5, kappa_R1: 0.7, kappa_R2: 1.3, q_tot: 1, K_Ao: 2, K_Ai: 2, q_Ao: 0.8, q_Ai: 0.3 }
    const expected = F.evaluate(limit.value.flux, values)
    const big = F.evaluate(full.value.flux, { ...values, kappa_R3: 1e9, kappa_R4: 1e9 })
    expect(big / expected).toBeCloseTo(1, 6)
  })

  it('reproduces the Kir channel’s flux with electrogenic steps (sodium homeostasis paper, eq 19)', () => {
    // States 1 = E_iK, 2 = E_oK, 3 = E; κ1 and κ3 (binding) fast; each reaction carries z_k F u/RT terms
    const zFu = (z: string) => times(ci(z), ci('F'), ci('u_m'))
    const equations = bondGraph([
      { name: 'R1', forward: { 3: 1, Ki: 1 }, reverse: { 1: 1 }, forwardPotential: zFu('z1'), reversePotential: zFu('z2') },
      { name: 'R2', forward: { 1: 1 }, reverse: { 2: 1 }, forwardPotential: zFu('z3'), reversePotential: zFu('z4') },
      { name: 'R3', forward: { 2: 1 }, reverse: { 3: 1, Ko: 1 }, forwardPotential: zFu('z5'), reversePotential: zFu('z6') },
    ])
    const result = solveCycle(equations, { cycleStates: ['q_1', 'q_2', 'q_3'], flux: 'v_R2', fast: ['kappa_R1', 'kappa_R3'], opaque: ['u_m'] })
    if (!result.ok) throw new Error(result.issues[0]!.message)
    const { flux, atoms } = result.value
    const values: Record<string, number> = {
      K_1: 2,
      K_2: 0.7,
      K_3: 1.4,
      K_Ki: 0.0239,
      K_Ko: 0.0239,
      q_Ki: 140,
      q_Ko: 4.5,
      kappa_R2: 0.9,
      q_tot: 1,
      z1: 0.1,
      z2: 0.05,
      z3: 0.2,
      z4: -0.3,
      z5: 0.15,
      z6: -0.4,
      F: 96485,
      RT: 2578.73,
      u_m: -0.06,
    }
    const e = (z: string) => Math.exp((values[z]! * values.F! * values.u_m!) / values.RT!)
    const qi = values.K_Ki! * values.q_Ki!
    const qo = values.K_Ko! * values.q_Ko!
    const expected =
      (values.q_tot! * values.K_3! * values.kappa_R2! * (qi * e('z1') / e('z2') * e('z3') - qo * e('z6') / e('z5') * e('z4'))) /
      (1 + (values.K_3! / values.K_1!) * qi * (e('z1') / e('z2')) + (values.K_3! / values.K_2!) * qo * (e('z6') / e('z5')))
    expect(F.evaluate(flux, atomValues(atoms, values))).toBeCloseTo(expected, 12)
  })

  it('reports fast rate constants that aren’t in the cycle', () => {
    const r = solveCycle(glut2(), { cycleStates: ['q_1', 'q_2', 'q_3', 'q_4'], flux: 'v_R2', fast: ['kappa_X'] })
    expect(r.ok ? undefined : r.issues[0]!.code).toBe('unknown-variable')
  })

  it('takes the fast limit of a single fast step', () => {
    // A two-state cycle with a fast step: q1 ⇌ q2 fast (with an exchange A), q2 → q1 slow
    const equations = bondGraph([
      { name: 'Rf', forward: { 1: 1, A: 1 }, reverse: { 2: 1 } },
      { name: 'Rs', forward: { 2: 1 }, reverse: { 1: 1, B: 1 } },
    ])
    const r = solveCycle(equations, { cycleStates: ['q_1', 'q_2'], flux: 'v_Rf', fast: ['kappa_Rf'] })
    if (!r.ok) throw new Error(r.issues[0]!.message)
    // At equilibrium K1 q1 K_A q_A = K2 q2; flux = κs(K2 q2 − K1 q1 K_B q_B)
    const values: Record<string, number> = { K_1: 2, K_2: 3, K_A: 1.5, K_B: 0.5, q_A: 0.4, q_B: 0.9, kappa_Rs: 1.1, q_tot: 1 }
    const q2PerQ1 = (values.K_1! * values.K_A! * values.q_A!) / values.K_2!
    const q1 = 1 / (1 + q2PerQ1)
    const expected = values.kappa_Rs! * (values.K_2! * q2PerQ1 * q1 - values.K_1! * q1 * values.K_B! * values.q_B!)
    expect(F.evaluate(r.value.flux, values)).toBeCloseTo(expected, 12)
  })
})
