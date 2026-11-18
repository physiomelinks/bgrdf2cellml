/**
 * Products of variables to integer powers, positive or negative, and how to
 * write them, with exp atoms combined into one `exp(…)`.
 */

import { cn, ci, divide, exp, minus, plus, power, times, type AstNode } from '@physiomelinks/math-ast'

import type { ExpAtom, Laurent } from './normalise.js'
import type { Monomial } from './poly/polynomial.js'
import * as Q from './poly/rational.js'
import type { Rational } from './poly/rational.js'

export type { Laurent }

const byName = (a: string, b: string): number => (a < b ? -1 : a > b ? 1 : 0)

export function laurent(powers: Iterable<readonly [string, number]>): Laurent {
  const merged = new Map<string, number>()
  for (const [v, e] of powers) merged.set(v, (merged.get(v) ?? 0) + e)
  return [...merged].filter(([, e]) => e !== 0).sort(([a], [b]) => byName(a, b))
}

export const fromMonomial = (m: Monomial): Laurent => m

export const multiply = (a: Laurent, b: Laurent): Laurent => laurent([...a, ...b])
export const inverse = (a: Laurent): Laurent => a.map(([v, e]) => [v, -e] as const)
export const divideLaurent = (a: Laurent, b: Laurent): Laurent => multiply(a, inverse(b))
export const pick = (a: Laurent, keep: (v: string) => boolean): Laurent => a.filter(([v]) => keep(v))
export const key = (a: Laurent): string => (a.length === 0 ? '1' : a.map(([v, e]) => `${v}^${e}`).join('*'))
export const isOne = (a: Laurent): boolean => a.length === 0

/** Only nonnegative powers: as a polynomial monomial. */
export function asMonomial(a: Laurent): Monomial | undefined {
  return a.every(([, e]) => e > 0) ? a : undefined
}

export function rationalAst(r: Rational): AstNode {
  return r.d === 1n ? cn(Number(r.n)) : divide(cn(Number(r.n)), cn(Number(r.d)))
}

function product(factors: AstNode[]): AstNode {
  return factors.length === 0 ? cn(1) : factors.length === 1 ? factors[0]! : times(...factors)
}

/** The exponent `Σ (e/L)·base` of a product of atoms, or `undefined` if there are none. */
export function atomExponent(m: Laurent, atoms: ReadonlyMap<string, ExpAtom>): AstNode | undefined {
  const terms: { r: Rational; base: AstNode }[] = []
  // Atoms with the same base would already be one atom, so each gives one term
  for (const [v, e] of m) {
    const atom = atoms.get(v)
    if (atom === undefined) continue
    terms.push({ r: Q.rational(e, atom.denominator), base: laurentAst(atom.base) })
  }
  if (terms.length === 0) return undefined
  let result: AstNode | undefined
  for (const { r, base } of terms) {
    const magnitude = Q.sign(r) < 0 ? Q.neg(r) : r
    const term = Q.isOne(magnitude) ? base : times(rationalAst(magnitude), base)
    if (result === undefined) result = Q.sign(r) < 0 ? minus(term) : term
    else if (Q.sign(r) < 0) result = minus(result, term)
    else result = result.type === 'Add' ? plus(...result.children, term) : plus(result, term)
  }
  return result
}

/**
 * A product as an expression tree: positive powers over negative ones, with
 * any atoms written as one `exp(…)` factor.
 */
export function laurentAst(m: Laurent, atoms: ReadonlyMap<string, ExpAtom> = new Map()): AstNode {
  const plain = m.filter(([v]) => !atoms.has(v))
  const factor = ([v, e]: readonly [string, number]): AstNode => (Math.abs(e) === 1 ? ci(v) : power(ci(v), cn(Math.abs(e))))
  const top = plain.filter(([, e]) => e > 0).map(factor)
  const argument = atomExponent(m, atoms)
  if (argument !== undefined) top.push(exp(argument))
  const bottom = plain.filter(([, e]) => e < 0).map(factor)
  return bottom.length === 0 ? product(top) : divide(product(top), product(bottom))
}
