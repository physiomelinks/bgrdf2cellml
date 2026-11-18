/**
 * Sparse multivariate polynomials with exact rational coefficients.
 *
 * Variables are names (CellML identifiers, or atoms the normaliser makes up).
 * A monomial is a product of variables to positive integer powers; a
 * polynomial maps each monomial to a nonzero coefficient. Values are
 * immutable: every operation returns a new polynomial.
 */

import * as Q from './rational.js'
import type { Rational } from './rational.js'

/** Variables and their positive powers, sorted by variable name. */
export type Monomial = readonly (readonly [string, number])[]

export interface Term {
  readonly monomial: Monomial
  readonly coefficient: Rational
}

const byName = (a: string, b: string): number => (a < b ? -1 : a > b ? 1 : 0)

/** A monomial's key: `a^2*b`, or `1` for the empty monomial. */
export function monomialKey(m: Monomial): string {
  return m.length === 0 ? '1' : m.map(([v, e]) => (e === 1 ? v : `${v}^${e}`)).join('*')
}

/** A monomial from a record or list of powers; zero powers are dropped. */
export function monomial(powers: Readonly<Record<string, number>> | Iterable<readonly [string, number]>): Monomial {
  const entries = Symbol.iterator in powers ? [...(powers as Iterable<readonly [string, number]>)] : Object.entries(powers)
  const merged = new Map<string, number>()
  for (const [v, e] of entries) {
    if (!Number.isInteger(e) || e < 0) throw new RangeError(`Not a nonnegative integer power: ${v}^${e}`)
    merged.set(v, (merged.get(v) ?? 0) + e)
  }
  return [...merged].filter(([, e]) => e !== 0).sort(([a], [b]) => byName(a, b))
}

export function multiplyMonomials(a: Monomial, b: Monomial): Monomial {
  const result: [string, number][] = []
  let i = 0
  let j = 0
  while (i < a.length || j < b.length) {
    const x = a[i]
    const y = b[j]
    if (y === undefined || (x !== undefined && x[0] < y[0])) {
      result.push([x![0], x![1]])
      i++
    } else if (x === undefined || y[0] < x[0]) {
      result.push([y[0], y[1]])
      j++
    } else {
      result.push([x[0], x[1] + y[1]])
      i++
      j++
    }
  }
  return result
}

/** `a/b` if `b` divides `a`, else `undefined`. */
export function divideMonomials(a: Monomial, b: Monomial): Monomial | undefined {
  const powers = new Map(a)
  for (const [v, e] of b) {
    const have = powers.get(v) ?? 0
    if (have < e) return undefined
    if (have === e) powers.delete(v)
    else powers.set(v, have - e)
  }
  return [...powers].sort(([x], [y]) => byName(x, y))
}

/** The greatest common divisor of two monomials: each variable to its lower power. */
export function gcdMonomials(a: Monomial, b: Monomial): Monomial {
  const other = new Map(b)
  const result: [string, number][] = []
  for (const [v, e] of a) {
    const f = other.get(v)
    if (f !== undefined) result.push([v, Math.min(e, f)])
  }
  return result
}

export const degreeOf = (m: Monomial, variable: string): number => m.find(([v]) => v === variable)?.[1] ?? 0

/** Lexicographic order on variable names: positive if `a` comes first (is greater). */
export function compareMonomials(a: Monomial, b: Monomial): number {
  let i = 0
  let j = 0
  for (;;) {
    const x = a[i]
    const y = b[j]
    if (x === undefined && y === undefined) return 0
    if (y === undefined) return 1
    if (x === undefined) return -1
    if (x[0] < y[0]) return 1
    if (y[0] < x[0]) return -1
    if (x[1] !== y[1]) return x[1] - y[1]
    i++
    j++
  }
}

export class Polynomial {
  /** Terms by monomial key; no zero coefficients. */
  private readonly map: ReadonlyMap<string, Term>

  private constructor(map: ReadonlyMap<string, Term>) {
    this.map = map
  }

  static readonly zero = new Polynomial(new Map())
  static readonly one = Polynomial.constant(Q.ONE)

  /** A sum of terms, collecting like ones. */
  static fromTerms(terms: Iterable<Term>): Polynomial {
    const map = new Map<string, Term>()
    for (const term of terms) {
      if (Q.isZero(term.coefficient)) continue
      const key = monomialKey(term.monomial)
      const existing = map.get(key)
      if (existing === undefined) map.set(key, term)
      else {
        const coefficient = Q.add(existing.coefficient, term.coefficient)
        if (Q.isZero(coefficient)) map.delete(key)
        else map.set(key, { monomial: existing.monomial, coefficient })
      }
    }
    return new Polynomial(map)
  }

  static constant(value: Rational | number): Polynomial {
    const coefficient = typeof value === 'number' ? Q.fromNumber(value) : value
    return Polynomial.fromTerms([{ monomial: [], coefficient }])
  }

  static variable(name: string, power = 1): Polynomial {
    return Polynomial.fromTerms([{ monomial: monomial([[name, power]]), coefficient: Q.ONE }])
  }

  static term(coefficient: Rational | number, m: Monomial): Polynomial {
    return Polynomial.fromTerms([{ monomial: m, coefficient: typeof coefficient === 'number' ? Q.fromNumber(coefficient) : coefficient }])
  }

  /** The terms, greatest monomial first (lexicographic on variable names). */
  terms(): Term[] {
    return [...this.map.values()].sort((a, b) => compareMonomials(b.monomial, a.monomial))
  }

  get size(): number {
    return this.map.size
  }

  isZero(): boolean {
    return this.map.size === 0
  }

  /** The constant's value if this is a constant (zero included), else `undefined`. */
  constantValue(): Rational | undefined {
    if (this.map.size === 0) return Q.ZERO
    if (this.map.size === 1) {
      const only = this.map.get('1')
      if (only !== undefined) return only.coefficient
    }
    return undefined
  }

  /** The coefficient of exactly this monomial (zero if absent). */
  coefficientOf(m: Monomial): Rational {
    return this.map.get(monomialKey(m))?.coefficient ?? Q.ZERO
  }

  /** Every variable that appears, sorted. */
  variables(): string[] {
    const names = new Set<string>()
    for (const { monomial: m } of this.map.values()) for (const [v] of m) names.add(v)
    return [...names].sort(byName)
  }

  has(variable: string): boolean {
    for (const { monomial: m } of this.map.values()) if (degreeOf(m, variable) > 0) return true
    return false
  }

  equals(other: Polynomial): boolean {
    if (this.map.size !== other.map.size) return false
    for (const [key, term] of this.map) {
      const o = other.map.get(key)
      if (o === undefined || !Q.equal(o.coefficient, term.coefficient)) return false
    }
    return true
  }

  add(other: Polynomial): Polynomial {
    if (other.isZero()) return this
    if (this.isZero()) return other
    return Polynomial.fromTerms([...this.map.values(), ...other.map.values()])
  }

  neg(): Polynomial {
    return Polynomial.fromTerms([...this.map.values()].map((t) => ({ monomial: t.monomial, coefficient: Q.neg(t.coefficient) })))
  }

  sub(other: Polynomial): Polynomial {
    return this.add(other.neg())
  }

  scale(factor: Rational | number): Polynomial {
    const f = typeof factor === 'number' ? Q.fromNumber(factor) : factor
    if (Q.isZero(f)) return Polynomial.zero
    if (Q.isOne(f)) return this
    return Polynomial.fromTerms([...this.map.values()].map((t) => ({ monomial: t.monomial, coefficient: Q.mul(t.coefficient, f) })))
  }

  /** Multiply by a monomial. */
  shift(m: Monomial): Polynomial {
    if (m.length === 0) return this
    return Polynomial.fromTerms([...this.map.values()].map((t) => ({ monomial: multiplyMonomials(t.monomial, m), coefficient: t.coefficient })))
  }

  mul(other: Polynomial): Polynomial {
    if (this.isZero() || other.isZero()) return Polynomial.zero
    const terms: Term[] = []
    for (const a of this.map.values())
      for (const b of other.map.values())
        terms.push({ monomial: multiplyMonomials(a.monomial, b.monomial), coefficient: Q.mul(a.coefficient, b.coefficient) })
    return Polynomial.fromTerms(terms)
  }

  pow(k: number): Polynomial {
    if (!Number.isInteger(k) || k < 0) throw new RangeError(`Not a nonnegative integer power: ${k}`)
    let result = Polynomial.one
    let base: Polynomial = this
    while (k > 0) {
      if (k & 1) result = result.mul(base)
      base = base.mul(base)
      k >>= 1
    }
    return result
  }

  /** The highest power of `variable`, or -1 for the zero polynomial. */
  degree(variable: string): number {
    let d = -1
    for (const { monomial: m } of this.map.values()) d = Math.max(d, degreeOf(m, variable))
    return d
  }

  /** The highest total power of the given variables over all terms, or -1 for zero. */
  degreeIn(variables: Iterable<string>): number {
    const names = new Set(variables)
    let d = -1
    for (const { monomial: m } of this.map.values()) d = Math.max(d, m.reduce((s, [v, e]) => (names.has(v) ? s + e : s), 0))
    return d
  }

  /**
   * Split by the powers of the given variables: each key is a monomial in
   * those variables, each value its coefficient, a polynomial in the rest.
   */
  coefficientsIn(variables: Iterable<string>): Map<string, { monomial: Monomial; coefficient: Polynomial }> {
    const names = new Set(variables)
    const groups = new Map<string, { monomial: Monomial; terms: Term[] }>()
    for (const t of this.map.values()) {
      const inside = t.monomial.filter(([v]) => names.has(v))
      const outside = t.monomial.filter(([v]) => !names.has(v))
      const key = monomialKey(inside)
      const group = groups.get(key) ?? { monomial: inside, terms: [] }
      group.terms.push({ monomial: outside, coefficient: t.coefficient })
      groups.set(key, group)
    }
    return new Map([...groups].map(([key, g]) => [key, { monomial: g.monomial, coefficient: Polynomial.fromTerms(g.terms) }]))
  }

  /** Replace each named variable by a polynomial. */
  substitute(values: ReadonlyMap<string, Polynomial> | Readonly<Record<string, Polynomial>>): Polynomial {
    const lookup = values instanceof Map ? values : new Map(Object.entries(values))
    if (lookup.size === 0) return this
    let result = Polynomial.zero
    const powers = new Map<string, Polynomial[]>()
    const power = (v: string, e: number): Polynomial => {
      const cache = powers.get(v) ?? [Polynomial.one]
      powers.set(v, cache)
      while (cache.length <= e) cache.push(cache[cache.length - 1]!.mul(lookup.get(v)!))
      return cache[e]!
    }
    for (const t of this.map.values()) {
      let product = Polynomial.term(t.coefficient, t.monomial.filter(([v]) => !lookup.has(v)))
      for (const [v, e] of t.monomial) if (lookup.has(v)) product = product.mul(power(v, e))
      result = result.add(product)
    }
    return result
  }

  /** The value at the given numbers; throws if a variable has none. */
  evaluate(values: ReadonlyMap<string, number> | Readonly<Record<string, number>>): number {
    const lookup = values instanceof Map ? values : new Map(Object.entries(values))
    let sum = 0
    for (const t of this.map.values()) {
      let product = Q.toNumber(t.coefficient)
      for (const [v, e] of t.monomial) {
        const x = lookup.get(v)
        if (x === undefined) throw new Error(`No value for ${v}`)
        product *= x ** e
      }
      sum += product
    }
    return sum
  }

  /** The greatest term in lexicographic order, or `undefined` for zero. */
  leadingTerm(): Term | undefined {
    let lead: Term | undefined
    for (const t of this.map.values()) if (lead === undefined || compareMonomials(t.monomial, lead.monomial) > 0) lead = t
    return lead
  }

  /**
   * The monomial and rational content: the largest monomial and positive
   * rational dividing every term (1 for zero).
   */
  content(): { monomial: Monomial; coefficient: Rational } {
    let m: Monomial | undefined
    let c: Rational = Q.ZERO
    for (const t of this.map.values()) {
      m = m === undefined ? t.monomial : gcdMonomials(m, t.monomial)
      c = Q.gcdRational(c, t.coefficient)
    }
    return { monomial: m ?? [], coefficient: Q.isZero(c) ? Q.ONE : c }
  }

  /** Divide by a monomial (and coefficient) that divides every term. */
  divideTerm(m: Monomial, coefficient: Rational = Q.ONE): Polynomial {
    return Polynomial.fromTerms(
      [...this.map.values()].map((t) => {
        const q = divideMonomials(t.monomial, m)
        if (q === undefined) throw new Error(`${monomialKey(m)} doesn't divide ${monomialKey(t.monomial)}`)
        return { monomial: q, coefficient: Q.div(t.coefficient, coefficient) }
      }),
    )
  }

  /**
   * `this/divisor` when it divides exactly, else `undefined`, by the
   * multivariate division algorithm in lexicographic order (exact division
   * always leaves no remainder in any monomial order).
   */
  divideExact(divisor: Polynomial): Polynomial | undefined {
    if (divisor.isZero()) throw new RangeError('Division by the zero polynomial')
    const lead = divisor.leadingTerm()!
    if (divisor.size === 1) {
      const terms: Term[] = []
      for (const t of this.map.values()) {
        const q = divideMonomials(t.monomial, lead.monomial)
        if (q === undefined) return undefined
        terms.push({ monomial: q, coefficient: Q.div(t.coefficient, lead.coefficient) })
      }
      return Polynomial.fromTerms(terms)
    }
    let remainder: Polynomial = this
    const quotient: Term[] = []
    while (!remainder.isZero()) {
      const r = remainder.leadingTerm()!
      const q = divideMonomials(r.monomial, lead.monomial)
      if (q === undefined) return undefined
      const term: Term = { monomial: q, coefficient: Q.div(r.coefficient, lead.coefficient) }
      quotient.push(term)
      remainder = remainder.sub(divisor.mul(Polynomial.fromTerms([term])))
    }
    return Polynomial.fromTerms(quotient)
  }

  toString(): string {
    if (this.isZero()) return '0'
    return this.terms()
      .map(({ monomial: m, coefficient: c }, i) => {
        const negative = Q.sign(c) < 0
        const magnitude = negative ? Q.neg(c) : c
        const body = m.length === 0 ? Q.toString(magnitude) : Q.isOne(magnitude) ? monomialKey(m) : `${Q.toString(magnitude)}*${monomialKey(m)}`
        return i === 0 ? (negative ? `-${body}` : body) : negative ? ` - ${body}` : ` + ${body}`
      })
      .join('')
  }
}
