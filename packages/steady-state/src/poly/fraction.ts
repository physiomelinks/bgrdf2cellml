/**
 * Quotients of polynomials. Only monomial factors are cancelled
 * automatically (and a denominator that divides exactly); there is no
 * general multivariate gcd, which the steady-state algebra never needs.
 */

import * as Q from './rational.js'
import type { Rational } from './rational.js'
import { Polynomial, gcdMonomials } from './polynomial.js'

export interface Fraction {
  readonly num: Polynomial
  readonly den: Polynomial
}

/** `num/den` with common monomial and rational factors cancelled, and a positive leading denominator coefficient. */
export function fraction(num: Polynomial, den: Polynomial = Polynomial.one): Fraction {
  if (den.isZero()) throw new RangeError('A fraction with a zero denominator')
  if (num.isZero()) return { num: Polynomial.zero, den: Polynomial.one }
  const constant = den.constantValue()
  if (constant !== undefined) return { num: num.scale(Q.div(Q.ONE, constant)), den: Polynomial.one }
  const exact = num.divideExact(den)
  if (exact !== undefined) return { num: exact, den: Polynomial.one }
  const a = num.content()
  const b = den.content()
  const m = gcdMonomials(a.monomial, b.monomial)
  // Scale so the denominator's content is 1 and its leading coefficient positive
  const lead = den.leadingTerm()!.coefficient
  const scale: Rational = Q.sign(lead) < 0 ? Q.neg(b.coefficient) : b.coefficient
  return { num: num.divideTerm(m, scale), den: den.divideTerm(m, scale) }
}

export const fromPolynomial = (p: Polynomial): Fraction => ({ num: p, den: Polynomial.one })

export function add(a: Fraction, b: Fraction): Fraction {
  if (a.den.equals(b.den)) return fraction(a.num.add(b.num), a.den)
  return fraction(a.num.mul(b.den).add(b.num.mul(a.den)), a.den.mul(b.den))
}

export const neg = (a: Fraction): Fraction => ({ num: a.num.neg(), den: a.den })
export const sub = (a: Fraction, b: Fraction): Fraction => add(a, neg(b))
export const mul = (a: Fraction, b: Fraction): Fraction => fraction(a.num.mul(b.num), a.den.mul(b.den))

export function div(a: Fraction, b: Fraction): Fraction {
  if (b.num.isZero()) throw new RangeError('Division by zero')
  return fraction(a.num.mul(b.den), a.den.mul(b.num))
}

export function pow(a: Fraction, k: number): Fraction {
  if (!Number.isInteger(k)) throw new RangeError(`Not an integer power: ${k}`)
  return k >= 0 ? fraction(a.num.pow(k), a.den.pow(k)) : fraction(a.den.pow(-k), a.num.pow(-k))
}

export const isPolynomial = (a: Fraction): boolean => a.den.constantValue() !== undefined

export function evaluate(a: Fraction, values: ReadonlyMap<string, number> | Readonly<Record<string, number>>): number {
  return a.num.evaluate(values) / a.den.evaluate(values)
}

export const toString = (a: Fraction): string => (isPolynomial(a) ? a.num.toString() : `(${a.num.toString()})/(${a.den.toString()})`)
