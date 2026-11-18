/**
 * Defaults for reducing a model: its conserved cycles, their held species,
 * which reactions to take as fast, and the flux to report.
 *
 * A cycle is a set of states whose rates of change sum to zero (a 0/1 vector
 * in the left null space of the stoichiometry) and which the cycle solver
 * accepts, so a glucose moiety (free plus bound glucose, also conserved)
 * isn't mistaken for one. Its fast reactions are suggested as the ones that
 * bind or release a held species, as long as one reaction stays slow.
 */

import { cn, ci, divide, times, walk, type AstNode, type Equation } from '@physiomelinks/math-ast'

import { solveCycle, splitEquations } from './cycle.js'
import type { SteadyStateSpec } from './derive.js'
import { SteadyStateError } from './issues.js'
import { dependents, expandLinear } from './linear.js'
import { normalise } from './normalise.js'
import * as Q from './poly/rational.js'
import type { Rational } from './poly/rational.js'

export interface SuggestOptions {
  /** Parameter values, used to recognise transformer ratios as multiples of F. */
  readonly values?: Readonly<Record<string, number>>
  /** Parameters that are transformer ratios (bg2cellml's manifest role `ratio`). */
  readonly ratios?: Iterable<string>
  /** The name of Faraday's constant in the model; `F`. */
  readonly faraday?: string
}

export interface SuggestedCycle {
  readonly spec: SteadyStateSpec
  /** The cycle's reactions, by flux variable. */
  readonly reactions: readonly string[]
  /** The held species. */
  readonly solutes: readonly string[]
}

/** Rational row reduction; returns a basis of the null space of `M` (columns are unknowns). */
function nullSpace(M: Rational[][], columns: number): Rational[][] {
  const m = M.map((r) => [...r])
  const pivots: number[] = []
  let row = 0
  for (let c = 0; c < columns && row < m.length; c++) {
    const p = m.findIndex((r, i) => i >= row && !Q.isZero(r[c]!))
    if (p < 0) continue
    ;[m[row], m[p]] = [m[p]!, m[row]!]
    const pivot = m[row]![c]!
    m[row] = m[row]!.map((x) => Q.div(x, pivot))
    for (let i = 0; i < m.length; i++)
      if (i !== row && !Q.isZero(m[i]![c]!)) {
        const f = m[i]![c]!
        m[i] = m[i]!.map((x, j) => Q.sub(x, Q.mul(f, m[row]![j]!)))
      }
    pivots.push(c)
    row++
  }
  const free = [...Array(columns).keys()].filter((c) => !pivots.includes(c))
  return free.map((f) => {
    const x = new Array<Rational>(columns).fill(Q.ZERO)
    x[f] = Q.ONE
    pivots.forEach((c, i) => (x[c] = Q.neg(m[i]![f]!)))
    return x
  })
}

/** A simple fraction a/b (b ≤ 12) within a relative 1e-5 of `x` (values are often rounded), if there is one. */
function simpleFraction(x: number): Rational | undefined {
  for (let b = 1; b <= 12; b++) {
    const a = Math.round(x * b)
    if (a !== 0 && Math.abs(a / b - x) <= 1e-5 * Math.abs(x)) return Q.rational(a, b)
  }
  return undefined
}

const kappaOf = (definition: AstNode | undefined, isParameter: (name: string) => boolean): string | undefined => {
  if (definition?.type === 'Divide') return kappaOf(definition.numerator, isParameter)
  if (definition?.type !== 'Multiply') return undefined
  const names = definition.children.flatMap((c) => (c.type === 'Identifier' && isParameter(c.name) ? [c.name] : []))
  return names.find((n) => n.startsWith('kappa')) ?? names[0]
}

export function suggestSpecs(equations: readonly Equation[], options: SuggestOptions = {}): SuggestedCycle[] {
  const { definitions, rates } = splitEquations(equations)
  const states = [...rates.keys()]
  const isDependent = dependents(definitions, states)
  const isParameter = (name: string) => !definitions.has(name) && !rates.has(name) && name !== 't'

  // Each state's rate as integer-ish multiples of reaction fluxes
  const rows = new Map<string, Map<string, Rational>>()
  for (const q of states) {
    const terms = expandLinear(rates.get(q)!, definitions, isDependent)
    const row = new Map<string, Rational>()
    let regular = true
    for (const t of terms) {
      if (t.flux === undefined || states.includes(t.flux)) {
        regular = false
        break
      }
      let c: Rational | undefined
      try {
        const f = normalise([t.coefficient], { definitions: new Map() }).fractions[0]!
        c = f.den.constantValue() === undefined ? undefined : f.num.constantValue()
        if (c !== undefined) c = Q.div(c, f.den.constantValue()!)
      } catch (e) {
        if (!(e instanceof SteadyStateError)) throw e
      }
      if (c === undefined) {
        regular = false
        break
      }
      row.set(t.flux, Q.add(row.get(t.flux) ?? Q.ZERO, c))
    }
    if (regular && row.size > 0) rows.set(q, row)
  }
  const regularStates = [...rows.keys()]
  const reactions = [...new Set([...rows.values()].flatMap((r) => [...r.keys()]))].sort()

  // 0/1 combinations of the left null space: states whose rates sum to zero
  const NT = reactions.map((r) => regularStates.map((q) => rows.get(q)!.get(r) ?? Q.ZERO))
  const basis = nullSpace(NT, regularStates.length)
  const candidates: number[][] = []
  const dimension = Math.min(basis.length, 6)
  const combos = 3 ** dimension
  for (let code = 1; code < combos; code++) {
    let k = code
    let x = new Array<Rational>(regularStates.length).fill(Q.ZERO)
    for (let b = 0; b < dimension; b++) {
      const coefficient = (k % 3) - 1
      k = Math.floor(k / 3)
      if (coefficient !== 0) x = x.map((v, i) => Q.add(v, Q.mul(Q.rational(coefficient), basis[b]![i]!)))
    }
    if (x.every((v) => Q.isZero(v) || Q.isOne(v)) && x.some((v) => Q.isOne(v))) candidates.push(x.flatMap((v, i) => (Q.isOne(v) ? [i] : [])))
  }
  candidates.sort((a, b) => a.length - b.length)

  const faraday = options.faraday ?? 'F'
  const F = options.values?.[faraday] ?? 96485.33212
  const relations: Record<string, AstNode> = {}
  for (const p of options.ratios ?? []) {
    const v = options.values?.[p]
    if (v === undefined) continue
    const times_ = simpleFraction(v / F)
    const over = simpleFraction(v * F)
    const r = (q: Rational) => (q.d === 1n ? cn(Number(q.n)) : divide(cn(Number(q.n)), cn(Number(q.d))))
    if (times_ !== undefined) relations[p] = times(r(times_), ci(faraday))
    else if (over !== undefined) relations[p] = divide(r(over), ci(faraday))
  }
  const withRelations = Object.keys(relations).length > 0 ? { relations } : {}

  const chosen: SuggestedCycle[] = []
  const used = new Set<string>()
  for (const candidate of candidates) {
    const cycleStates = candidate.map((i) => regularStates[i]!)
    if (cycleStates.length < 2 || cycleStates.some((q) => used.has(q))) continue
    const cycleReactions = reactions.filter((r) => cycleStates.some((q) => !Q.isZero(rows.get(q)!.get(r) ?? Q.ZERO)))
    // A cycle closes on itself; a moiety (free plus bound glucose, say) is a chain with fewer reactions than states
    if (cycleReactions.length < cycleStates.length) continue
    const solutes = states.filter((q) => !cycleStates.includes(q) && cycleReactions.some((r) => rows.get(q) !== undefined ? !Q.isZero(rows.get(q)!.get(r) ?? Q.ZERO) : false))
    const binding = cycleReactions.filter((r) => solutes.some((q) => !Q.isZero(rows.get(q)?.get(r) ?? Q.ZERO)))
    const slow = cycleReactions.filter((r) => !binding.includes(r))
    const kappas = binding.map((r) => kappaOf(definitions.get(r), isParameter))
    const fast = slow.length > 0 && kappas.every((k) => k !== undefined) ? (kappas as string[]) : []
    const flux = (fast.length > 0 ? slow : cycleReactions)[0]!

    // Membrane potentials and similar: defined from a state outside the cycle, with no logarithm
    const opaque = [...definitions].flatMap(([name, d]) => {
      const nodes = [...walk(d)]
      const fromState = nodes.some((n) => n.type === 'Identifier' && rates.has(n.name) && !cycleStates.includes(n.name))
      const plain = nodes.every((n) => n.type !== 'FunctionCall')
      const simple = nodes.every((n) => n.type !== 'Identifier' || !definitions.has(n.name))
      return fromState && plain && simple ? [name] : []
    })

    const base: SteadyStateSpec = { cycleStates, flux, ...(opaque.length > 0 ? { opaque } : {}), ...withRelations }
    for (const spec of fast.length > 0 ? [{ ...base, fast }, base] : [base]) {
      if (solveCycle([...equations, ...Object.entries(relations).map(([n, v]) => ({ type: 'Equal' as const, left: ci(n), right: v }))], spec).ok) {
        chosen.push({ spec, reactions: cycleReactions, solutes })
        for (const q of cycleStates) used.add(q)
        break
      }
    }
  }
  return chosen
}
