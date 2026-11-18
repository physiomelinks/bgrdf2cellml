/**
 * A numeric check of a derived flux, independent of the algebra: at random
 * values, solve the model's own cycle equations for its steady state with
 * plain floating-point linear algebra and compare the flux.
 */

import { ci, eq, walk, type AstNode, type Equation } from '@physiomelinks/math-ast'

import { splitEquations } from './cycle.js'
import type { SteadyStateResult, SteadyStateSpec } from './derive.js'
import { evaluateAst } from './evaluate.js'
import { atomValues } from './normalise.js'
import * as F from './poly/fraction.js'

export interface VerifyOptions {
  readonly samples?: number
  /** Values to use rather than random ones. */
  readonly values?: Readonly<Record<string, number>>
  /** The value fast rate constants get, relative to the others; 1e7 if not given. */
  readonly fastScale?: number
  /** A seed, for repeatable samples. */
  readonly seed?: number
}

export interface VerifyReport {
  /** The largest relative difference between the derived flux and the numeric one. */
  readonly maxRelativeError: number
  /** The same for the common form's expression, when there is one. */
  readonly commonFormError?: number
  readonly samples: number
}

function solveLinear(A: number[][], b: number[]): number[] {
  const n = b.length
  const m = A.map((row, i) => [...row, b[i]!])
  for (let k = 0; k < n; k++) {
    let p = k
    for (let i = k + 1; i < n; i++) if (Math.abs(m[i]![k]!) > Math.abs(m[p]![k]!)) p = i
    ;[m[k], m[p]] = [m[p]!, m[k]!]
    for (let i = k + 1; i < n; i++) {
      const f = m[i]![k]! / m[k]![k]!
      for (let j = k; j <= n; j++) m[i]![j]! -= f * m[k]![j]!
    }
  }
  const x = new Array<number>(n).fill(0)
  for (let i = n - 1; i >= 0; i--) {
    let s = m[i]![n]!
    for (let j = i + 1; j < n; j++) s -= m[i]![j]! * x[j]!
    x[i] = s / m[i]![i]!
  }
  return x
}

/** A small seeded generator (mulberry32). */
function random(seed: number): () => number {
  let a = seed >>> 0
  return () => {
    a = (a + 0x6d2b79f5) >>> 0
    let t = a
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}

export function verify(equations: readonly Equation[], spec: SteadyStateSpec, result: SteadyStateResult, options: VerifyOptions = {}): VerifyReport {
  if (result.solution === undefined) throw new Error('Nothing to verify: the cycle wasn’t solved')
  const all = [...equations, ...Object.entries(spec.relations ?? {}).map(([name, value]) => eq(ci(name), value))]
  const { definitions, rates } = splitEquations(all)
  const states = [...spec.cycleStates]
  const total = result.solution.total
  const fast = new Set(spec.fast ?? [])
  const fluxNode: AstNode = typeof spec.flux === 'string' ? ci(spec.flux) : spec.flux

  // The free variables: neither defined nor a cycle state
  const free = new Set<string>()
  for (const node of [...definitions.values(), ...states.map((q) => rates.get(q)!), fluxNode])
    for (const n of walk(node)) if (n.type === 'Identifier' && !definitions.has(n.name) && !states.includes(n.name)) free.add(n.name)
  if (result.common !== undefined) {
    const lumped = new Set([...result.common.variables, ...result.common.definitions].map(({ name }) => name))
    for (const { value } of [...result.common.variables, ...result.common.definitions])
      for (const n of walk(value)) if (n.type === 'Identifier' && !lumped.has(n.name) && !definitions.has(n.name)) free.add(n.name)
  }
  free.delete('t')

  const next = random(options.seed ?? 1)
  const samples = options.samples ?? 5
  let maxRelativeError = 0
  let commonFormError: number | undefined
  for (let s = 0; s < samples; s++) {
    const values: Record<string, number> = { [total]: 1 }
    for (const v of free) values[v] = options.values?.[v] ?? (fast.has(v) ? (options.fastScale ?? 1e7) : 0.5 + 1.5 * next())
    // Each column of A: the states' rates of change with one state 1 and the rest 0
    const A = states.map(() => new Array<number>(states.length).fill(0))
    states.forEach((_qi, i) => {
      const at = { ...values, ...Object.fromEntries(states.map((q, k) => [q, k === i ? 1 : 0])) }
      states.forEach((qj, j) => (A[j]![i] = evaluateAst(rates.get(qj)!, at, definitions)))
    })
    A[states.length - 1] = states.map(() => 1)
    const b = states.map((_q, j) => (j === states.length - 1 ? values[total]! : 0))
    const q = solveLinear(A, b)
    const numeric = evaluateAst(fluxNode, { ...values, ...Object.fromEntries(states.map((name, k) => [name, q[k]!])) }, definitions)
    // Opaque variables the derived flux keeps, such as a membrane potential, from their definitions
    const env: Record<string, number> = { ...values }
    const sol = result.solution
    for (const v of new Set([...sol.flux.num.variables(), ...sol.flux.den.variables(), ...[...sol.atoms.values()].flatMap((a) => a.base.map(([x]) => x))]))
      if (env[v] === undefined && definitions.has(v)) env[v] = evaluateAst(ci(v), { ...values, ...Object.fromEntries(states.map((name, k) => [name, q[k]!])) }, definitions)
    const derived = F.evaluate(sol.flux, atomValues(sol.atoms, env))
    const scale = Math.max(Math.abs(numeric), Math.abs(derived), Number.MIN_VALUE)
    maxRelativeError = Math.max(maxRelativeError, Math.abs(numeric - derived) / scale)
    if (result.common !== undefined) {
      const lumped = new Map<string, AstNode>([...definitions, ...[...result.common.variables, ...result.common.definitions].map(({ name, value }) => [name, value] as const)])
      const common = evaluateAst(result.common.flux, env, lumped)
      commonFormError = Math.max(commonFormError ?? 0, Math.abs(common - derived) / Math.max(Math.abs(derived), Number.MIN_VALUE))
    }
  }
  return { maxRelativeError, ...(commonFormError === undefined ? {} : { commonFormError }), samples }
}
