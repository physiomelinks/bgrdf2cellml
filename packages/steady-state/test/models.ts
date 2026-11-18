/**
 * Bond-graph models written the way bg2cellml writes them, for tests: each
 * species `u = RT ln(K q)`, each reaction `v = κ (exp(u_f/RT) − exp(u_r/RT))`
 * with `u_f` and `u_r` the stoichiometry-weighted sums of potentials, and
 * `dq/dt` the signed sum of the reactions' fluxes.
 */

import { ci, cn, divide, eq, exp, fn, minus, plus, times, type AstNode, type Equation } from '@physiomelinks/math-ast'

export interface Reaction {
  readonly name: string
  /** Species and their stoichiometry on each side. */
  readonly forward: Readonly<Record<string, number>>
  readonly reverse: Readonly<Record<string, number>>
  /** Extra potential added to each side, such as `z F u_m` (J/mol). */
  readonly forwardPotential?: AstNode
  readonly reversePotential?: AstNode
}

const derivative = (q: string): AstNode => ({ type: 'Derivative', expression: ci(q), variable: ci('t') })

const weighted = (side: Readonly<Record<string, number>>, extra?: AstNode): AstNode => {
  const terms: AstNode[] = Object.entries(side).map(([s, n]) => (n === 1 ? ci(`u_${s}`) : times(cn(n), ci(`u_${s}`))))
  if (extra !== undefined) terms.push(extra)
  return terms.length === 1 ? terms[0]! : plus(...terms)
}

/** The equations of a bond-graph model of these reactions, every species with `K_<s>` and every reaction `kappa_<r>`. */
export function bondGraph(reactions: readonly Reaction[]): Equation[] {
  const species = new Set<string>()
  for (const r of reactions) for (const s of [...Object.keys(r.forward), ...Object.keys(r.reverse)]) species.add(s)
  const equations: Equation[] = []
  for (const s of species) equations.push(eq(ci(`u_${s}`), times(ci('RT'), fn('ln', times(ci(`K_${s}`), ci(`q_${s}`))))))
  for (const r of reactions) {
    equations.push(eq(ci(`u_f_${r.name}`), weighted(r.forward, r.forwardPotential)))
    equations.push(eq(ci(`u_r_${r.name}`), weighted(r.reverse, r.reversePotential)))
    equations.push(
      eq(
        ci(`v_${r.name}`),
        times(ci(`kappa_${r.name}`), minus(exp(divide(ci(`u_f_${r.name}`), ci('RT'))), exp(divide(ci(`u_r_${r.name}`), ci('RT'))))),
      ),
    )
  }
  for (const s of species) {
    const terms: AstNode[] = []
    for (const r of reactions) {
      const n = (r.reverse[s] ?? 0) - (r.forward[s] ?? 0)
      if (n === 0) continue
      const v = ci(`v_${r.name}`)
      const scaled = Math.abs(n) === 1 ? v : times(cn(Math.abs(n)), v)
      terms.push(n > 0 ? scaled : minus(scaled))
    }
    equations.push(eq(derivative(`q_${s}`), terms.length === 1 ? terms[0]! : plus(...terms)))
  }
  return equations
}

/** Michaelis–Menten: E + S ⇌ C ⇌ E + P (Hunter et al. 2025, Fig 5). */
export const enzyme = (): Equation[] =>
  bondGraph([
    { name: 'R1', forward: { S: 1, E: 1 }, reverse: { C: 1 } },
    { name: 'R2', forward: { C: 1 }, reverse: { P: 1, E: 1 } },
  ])

/** GLUT2's alternating-access cycle (Hunter et al. 2025, Fig 6b), states 1–4 with glucose A outside (o) and inside (i). */
export const glut2 = (): Equation[] =>
  bondGraph([
    { name: 'R1', forward: { 4: 1 }, reverse: { 1: 1 } },
    { name: 'R2', forward: { 2: 1 }, reverse: { 3: 1 } },
    { name: 'R3', forward: { Ao: 1, 1: 1 }, reverse: { 2: 1 } },
    { name: 'R4', forward: { 3: 1 }, reverse: { 4: 1, Ai: 1 } },
  ])

/** Solve `A x = b` numerically by Gaussian elimination with partial pivoting. */
export function solveLinear(A: number[][], b: number[]): number[] {
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
