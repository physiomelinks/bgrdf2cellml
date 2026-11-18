/**
 * The steady state of a linear cycle: an enzyme or transporter whose states'
 * rates of change are linear in those states, with the other species held
 * constant.
 *
 * Setting the states' rates of change to zero gives `A·q = 0`, where each
 * entry of A is a polynomial in the parameters, the held species and the exp
 * atoms. The states sum to a constant total, so A is singular and one of its
 * rows is redundant; replacing it by the total gives a system Cramer's rule
 * solves exactly. Each state's share is then a cofactor of A along the
 * dropped row: the King–Altman (matrix-tree) weights of the cycle's diagram.
 */

import { ci, type AstNode, type Equation } from '@physiomelinks/math-ast'

import { fail, SteadyStateError, type Issue } from './issues.js'
import { normalise, type ExpAtom } from './normalise.js'
import * as F from './poly/fraction.js'
import type { Fraction } from './poly/fraction.js'
import { cofactor } from './poly/matrix.js'
import { Polynomial, gcdMonomials, monomialKey, type Monomial } from './poly/polynomial.js'

export interface CycleSpec {
  /** The cycle's states, whose rates of change are set to zero. */
  readonly cycleStates: readonly string[]
  /** The flux to report: a variable's name or an expression. */
  readonly flux: string | AstNode
  /** The name for the states' total amount; `q_tot` if not given. */
  readonly total?: string
  /** Variables to keep as they are rather than inline, such as a membrane potential. */
  readonly opaque?: readonly string[]
}

export interface CycleSolution {
  /** The flux at steady state, as a quotient of polynomials in parameters, held species, `total` and exp atoms. */
  readonly flux: Fraction
  /** Each state's weight: at steady state `q_j = total · weights[j] / Σ weights`. */
  readonly weights: ReadonlyMap<string, Polynomial>
  /** `A[j][i]`: the coefficient of state i in the (cleared) rate of change of state j. */
  readonly matrix: readonly (readonly Polynomial[])[]
  readonly total: string
  readonly atoms: ReadonlyMap<string, ExpAtom>
}

export type Result<T> = { readonly ok: true; readonly value: T } | { readonly ok: false; readonly issues: readonly Issue[] }

/** A model's equations split into definitions (`x = …`) and rates of change (`d(q)/dt = …`). */
export interface SplitEquations {
  readonly definitions: Map<string, AstNode>
  readonly rates: Map<string, AstNode>
}

export function splitEquations(equations: readonly Equation[]): SplitEquations {
  const definitions = new Map<string, AstNode>()
  const rates = new Map<string, AstNode>()
  for (const { left, right } of equations) {
    if (left.type === 'Identifier') definitions.set(left.name, right)
    else if (left.type === 'Derivative' && left.expression.type === 'Identifier') rates.set(left.expression.name, right)
  }
  return { definitions, rates }
}

/** The steady-state flux of a linear cycle, or the issues that stop it. */
export function solveCycle(equations: readonly Equation[], spec: CycleSpec): Result<CycleSolution> {
  try {
    return { ok: true, value: solve(equations, spec) }
  } catch (e) {
    if (e instanceof SteadyStateError) return { ok: false, issues: [e.issue] }
    throw e
  }
}

function solve(equations: readonly Equation[], spec: CycleSpec): CycleSolution {
  const { definitions, rates } = splitEquations(equations)
  const states = [...spec.cycleStates]
  const total = spec.total ?? 'q_tot'
  if (states.length < 2) fail('disconnected-cycle', 'A cycle needs at least two states')
  for (const q of states) if (!rates.has(q)) fail('unknown-variable', `${q} has no rate of change in the model`, q)
  const fluxNode = typeof spec.flux === 'string' ? ci(spec.flux) : spec.flux

  const { fractions, atoms } = normalise([...states.map((q) => rates.get(q)!), fluxNode], {
    definitions,
    states,
    ...(spec.opaque === undefined ? {} : { opaque: spec.opaque }),
  })
  const stateRates = fractions.slice(0, states.length)
  const flux = fractions[states.length]!

  // The states must be conserved: their rates of change sum to zero
  if (!stateRates.reduce(F.add).num.isZero())
    fail('no-conservation', `The rates of change of ${states.join(', ')} don't sum to zero, so their total isn't conserved`)

  // Each rate's numerator must be linear in the states; denominators mustn't involve them
  const linear = (f: Fraction, what: string): Map<string, Polynomial> => {
    for (const q of states)
      if (f.den.has(q)) fail('nonlinear-in-states', `${what} has ${q} in a denominator, so it isn't linear in the cycle's states`, q)
    const groups = f.num.coefficientsIn(states)
    const coefficients = new Map<string, Polynomial>()
    for (const { monomial: m, coefficient } of groups.values()) {
      const degree = m.reduce((s, [, e]) => s + e, 0)
      if (degree === 0) coefficients.set('', coefficient)
      else if (degree === 1) coefficients.set(m[0]![0], coefficient)
      else fail('nonlinear-in-states', `${what} has a term in ${monomialKey(m)}, so it isn't linear in the cycle's states`)
    }
    return coefficients
  }

  const matrix = stateRates.map((f, j) => {
    const coefficients = linear(f, `d${states[j]}/dt`)
    if (coefficients.has('')) fail('nonlinear-in-states', `d${states[j]}/dt has a term without any of the cycle's states`, states[j])
    return states.map((q) => coefficients.get(q) ?? Polynomial.zero)
  })

  // Cofactors along the last row: a state's weight, up to a common factor
  const last = states.length - 1
  let weights = states.map((_q, i) => cofactor(matrix, last, i))
  const sum = weights.reduce((a, b) => a.add(b))
  if (sum.isZero()) fail('disconnected-cycle', `The states ${states.join(', ')} don't form one connected cycle, so their steady state isn't unique`)

  // Cancel what every weight shares, and make the sum's leading coefficient positive
  let common: Monomial | undefined
  for (const w of weights) if (!w.isZero()) common = common === undefined ? w.content().monomial : gcdMonomials(common, w.content().monomial)
  const sign = sum.leadingTerm()!.coefficient.n < 0n ? -1 : 1
  weights = weights.map((w) => (w.isZero() ? w : w.divideTerm(common ?? [], { n: BigInt(sign), d: 1n })))
  const weightSum = weights.reduce((a, b) => a.add(b))

  // The flux, with q_j = total · w_j / Σw
  const fluxCoefficients = linear(flux, 'The flux')
  const q = Polynomial.variable(total)
  let numerator = (fluxCoefficients.get('') ?? Polynomial.zero).mul(weightSum)
  states.forEach((s, i) => {
    const c = fluxCoefficients.get(s)
    if (c !== undefined) numerator = numerator.add(c.mul(weights[i]!).mul(q))
  })
  return {
    flux: F.fraction(numerator, flux.den.mul(weightSum)),
    weights: new Map(states.map((s, i) => [s, weights[i]!])),
    matrix,
    total,
    atoms,
  }
}
