import { describe, expect, it } from 'vitest'

import { Polynomial, Q, cofactor, determinant, monomial } from '../src/index.js'

const x = Polynomial.variable('x')
const y = Polynomial.variable('y')
const z = Polynomial.variable('z')
const c = (v: number) => Polynomial.constant(v)

describe('Rational', () => {
  it('keeps lowest terms with a positive denominator', () => {
    expect(Q.rational(6, -4)).toEqual({ n: -3n, d: 2n })
    expect(Q.rational(0, 5)).toEqual(Q.ZERO)
  })

  it('reads decimals exactly', () => {
    expect(Q.fromDecimal('0.1')).toEqual(Q.rational(1, 10))
    expect(Q.fromDecimal('2.5e-3')).toEqual(Q.rational(1, 400))
    expect(Q.fromDecimal('-1.2E2')).toEqual(Q.rational(-120))
    expect(Q.fromNumber(0.726)).toEqual(Q.rational(363, 500))
    expect(() => Q.fromDecimal('abc')).toThrow(SyntaxError)
  })

  it('does arithmetic', () => {
    const a = Q.rational(1, 3)
    const b = Q.rational(1, 6)
    expect(Q.add(a, b)).toEqual(Q.rational(1, 2))
    expect(Q.sub(b, a)).toEqual(Q.rational(-1, 6))
    expect(Q.mul(a, b)).toEqual(Q.rational(1, 18))
    expect(Q.div(a, b)).toEqual(Q.rational(2))
    expect(Q.pow(Q.rational(2, 3), -2)).toEqual(Q.rational(9, 4))
    expect(Q.gcdRational(Q.rational(4, 3), Q.rational(6, 5))).toEqual(Q.rational(2, 15))
    expect(() => Q.div(a, Q.ZERO)).toThrow(RangeError)
  })
})

describe('Polynomial', () => {
  it('collects like terms and drops zeros', () => {
    expect(x.add(y).sub(x).equals(y)).toBe(true)
    expect(x.sub(x).isZero()).toBe(true)
    expect(x.add(y).pow(2).toString()).toBe('x^2 + 2*x*y + y^2')
  })

  it('multiplies and expands', () => {
    const p = x.add(c(1)).mul(x.sub(c(1)))
    expect(p.equals(x.pow(2).sub(c(1)))).toBe(true)
    expect(x.mul(y).mul(x).toString()).toBe('x^2*y')
  })

  it('knows degrees and coefficients', () => {
    const p = x.pow(2).mul(y).add(x.mul(z)).add(c(3))
    expect(p.degree('x')).toBe(2)
    expect(p.degreeIn(['x', 'y'])).toBe(3)
    const byX = p.coefficientsIn(['x'])
    expect(byX.get('x^2')!.coefficient.equals(y)).toBe(true)
    expect(byX.get('x')!.coefficient.equals(z)).toBe(true)
    expect(byX.get('1')!.coefficient.equals(c(3))).toBe(true)
    expect(p.variables()).toEqual(['x', 'y', 'z'])
  })

  it('substitutes and evaluates', () => {
    const p = x.pow(2).add(y)
    expect(p.substitute({ x: y.add(c(1)) }).equals(y.pow(2).add(y.scale(3)).add(c(1)))).toBe(true)
    expect(p.evaluate({ x: 3, y: 0.5 })).toBeCloseTo(9.5)
    expect(() => p.evaluate({ x: 1 })).toThrow(/No value for y/)
  })

  it('finds the content', () => {
    const p = x.pow(2).mul(y).scale(Q.rational(4, 3)).add(x.mul(y.pow(3)).scale(2))
    const { monomial: m, coefficient } = p.content()
    expect(m).toEqual(monomial({ x: 1, y: 1 }))
    expect(coefficient).toEqual(Q.rational(2, 3))
    expect(p.divideTerm(m, coefficient).equals(x.scale(2).add(y.pow(2).scale(3)))).toBe(true)
  })

  it('divides exactly, or says it cannot', () => {
    const a = x.add(y).mul(x.sub(z)).mul(y.add(c(2)))
    expect(a.divideExact(x.sub(z))!.equals(x.add(y).mul(y.add(c(2))))).toBe(true)
    expect(a.divideExact(x.mul(y))).toBeUndefined()
    expect(x.add(c(1)).divideExact(x.sub(c(1)))).toBeUndefined()
  })
})

describe('determinant', () => {
  it('matches the cofactor expansion of a symbolic 3×3', () => {
    const names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'].map((n) => Polynomial.variable(n))
    const [a, b, cc, d, e, f, g, h, i] = names as [Polynomial, ...Polynomial[]] & Polynomial[]
    const m = [
      [a, b, cc],
      [d, e, f],
      [g, h, i],
    ] as Polynomial[][]
    const expected = a!
      .mul(e!.mul(i!).sub(f!.mul(h!)))
      .sub(b!.mul(d!.mul(i!).sub(f!.mul(g!))))
      .add(cc!.mul(d!.mul(h!).sub(e!.mul(g!))))
    expect(determinant(m).equals(expected)).toBe(true)
  })

  it('pivots past zeros and handles singular matrices', () => {
    const zero = Polynomial.zero
    expect(
      determinant([
        [zero, x, c(0)],
        [y, c(0), c(0)],
        [c(0), c(0), z],
      ]).equals(x.mul(y).mul(z).neg()),
    ).toBe(true)
    expect(
      determinant([
        [x, y],
        [x.scale(2), y.scale(2)],
      ]).isZero(),
    ).toBe(true)
  })

  it('gives the 4-state cycle steady state by its cofactors (King–Altman)', () => {
    // A cycle 1→2→3→4→1 with forward rates k_i and reverse rates r_i
    const k = [1, 2, 3, 4].map((n) => Polynomial.variable(`k${n}`))
    const r = [1, 2, 3, 4].map((n) => Polynomial.variable(`r${n}`))
    const n = 4
    // A[j][i]: rate of change of q_j per unit q_i
    const A: Polynomial[][] = Array.from({ length: n }, () => Array.from({ length: n }, () => Polynomial.zero))
    for (let s = 0; s < n; s++) {
      const t = (s + 1) % n
      A[t]![s] = A[t]![s]!.add(k[s]!)
      A[s]![s] = A[s]![s]!.sub(k[s]!)
      A[s]![t] = A[s]![t]!.add(r[s]!)
      A[t]![t] = A[t]![t]!.sub(r[s]!)
    }
    const D = Array.from({ length: n }, (_, j) => cofactor(A, j, j))
    // Each state's weight is a sum of 4 spanning trees of 3 rates each
    for (const d of D) {
      expect(d.size).toBe(4)
      expect(d.degreeIn([...k, ...r].flatMap((p) => p.variables()))).toBe(3)
    }
    // And A·D = 0 at random values
    const values: Record<string, number> = {}
    for (const p of [...k, ...r]) values[p.variables()[0]!] = 0.5 + Math.random()
    for (let j = 0; j < n; j++) {
      const row = A[j]!.reduce((s, a, i) => s + a.evaluate(values) * D[i]!.evaluate(values), 0)
      expect(Math.abs(row)).toBeLessThan(1e-9)
    }
  })
})
