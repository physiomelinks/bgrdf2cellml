/**
 * Rapid equilibrium (Cha's method): the limit as some reactions' rate
 * constants κ go to infinity, such as a transporter's binding and unbinding.
 *
 * States joined by fast reactions equilibrate among themselves at once, so
 * each set of them is a pool whose members are in fixed proportions: those
 * that make every fast reaction's forward and reverse terms equal (the paper's
 * eqs 15–16 and 30–33). The slow reactions then move amount between pools,
 * and the pooled cycle is solved like any other. A flux through a fast
 * reaction is finite though its κ isn't: it is rewritten first, using the
 * states' steady-state balances, in terms of slow reactions only.
 */

import { fail } from './issues.js'
import * as F from './poly/fraction.js'
import type { Fraction } from './poly/fraction.js'
import { cofactor } from './poly/matrix.js'
import { Polynomial, gcdMonomials, type Monomial, type Term } from './poly/polynomial.js'
import * as Q from './poly/rational.js'
import type { Rational } from './poly/rational.js'

/** Split a polynomial into its terms with a fast rate constant and those without. */
function splitFast(p: Polynomial, fast: ReadonlySet<string>): { fast: Polynomial; slow: Polynomial } {
  const fastTerms: Term[] = []
  const slowTerms: Term[] = []
  for (const t of p.terms()) {
    const powers = t.monomial.filter(([v]) => fast.has(v))
    if (powers.length > 1 || (powers[0] !== undefined && powers[0][1] > 1))
      fail('unsupported-construct', 'A term with more than one fast rate constant, so the fast limit isn’t defined')
    ;(powers.length === 0 ? slowTerms : fastTerms).push(t)
  }
  return { fast: Polynomial.fromTerms(fastTerms), slow: Polynomial.fromTerms(slowTerms) }
}

/** Solve `M x = b` over the rationals, with free unknowns 0; `undefined` if inconsistent. */
function solveRational(M: Rational[][], b: Rational[]): Rational[] | undefined {
  const rows = M.length
  const cols = M[0]?.length ?? 0
  const m = M.map((row, i) => [...row, b[i]!])
  const pivots: number[] = []
  let r = 0
  for (let c = 0; c < cols && r < rows; c++) {
    const p = m.findIndex((row, i) => i >= r && !Q.isZero(row[c]!))
    if (p < 0) continue
    ;[m[r], m[p]] = [m[p]!, m[r]!]
    const pivot = m[r]![c]!
    m[r] = m[r]!.map((x) => Q.div(x, pivot))
    for (let i = 0; i < rows; i++) {
      if (i === r || Q.isZero(m[i]![c]!)) continue
      const f = m[i]![c]!
      m[i] = m[i]!.map((x, j) => Q.sub(x, Q.mul(f, m[r]![j]!)))
    }
    pivots.push(c)
    r++
  }
  for (let i = r; i < rows; i++) if (!Q.isZero(m[i]![cols]!)) return undefined
  const x = new Array<Rational>(cols).fill(Q.ZERO)
  pivots.forEach((c, i) => (x[c] = m[i]![cols]!))
  return x
}

/**
 * A flux equal, at steady state, to `flux` but without any fast rate
 * constant: `flux − Σ α_k dq_k/dt` with constant α chosen to cancel every
 * fast reaction's terms. `rows` are the states' rates of change.
 */
export function withoutFastReactions(flux: Polynomial, rows: readonly Polynomial[], fast: ReadonlySet<string>): Polynomial {
  const kappas = [...fast].filter((k) => flux.has(k) || rows.some((r) => r.has(k)))
  if (!kappas.some((k) => flux.has(k))) return flux
  // For each fast reaction r: its rate as it appears in the rows, and each row's and the flux's multiple of it
  const M: Rational[][] = []
  const b: Rational[] = []
  for (const kappa of kappas) {
    const part = (p: Polynomial): Polynomial => p.coefficientsIn([kappa]).get(kappa)?.coefficient ?? Polynomial.zero
    const rowParts = rows.map(part)
    const unit = rowParts.find((p) => !p.isZero())
    const fluxPart = part(flux)
    if (unit === undefined) {
      if (!fluxPart.isZero()) fail('unsupported-construct', `The flux uses the fast rate constant ${kappa}, which no cycle state's balance has`)
      continue
    }
    const multiple = (p: Polynomial): Rational => {
      if (p.isZero()) return Q.ZERO
      const q = p.divideExact(unit)?.constantValue()
      if (q === undefined) return fail('unsupported-construct', `The reaction with ${kappa} doesn't enter the states' balances as one rate`)
      return q
    }
    M.push(rowParts.map(multiple))
    b.push(multiple(fluxPart))
  }
  // α_k for each row k, with Σ_k α_k · s_kr = β_r for each fast reaction r
  const alpha = M.length === 0 ? [] : solveRational(M, b)
  if (alpha === undefined) return fail('not-thermodynamically-consistent', 'The flux through the fast reactions can’t be found from the slow ones')
  let result = flux
  alpha.forEach((a, k) => {
    if (!Q.isZero(a)) result = result.sub(rows[k]!.scale(a))
  })
  for (const kappa of kappas) if (result.has(kappa)) fail('unsupported-construct', `The flux still depends on the fast rate constant ${kappa}`)
  return result
}

/**
 * The states' steady-state weights in the fast limit: `q_j = total · w_j / Σ w`.
 * `matrix[j][i]` is the coefficient of state i in state j's rate of change.
 */
export function rapidEquilibriumWeights(matrix: readonly (readonly Polynomial[])[], fast: ReadonlySet<string>): Polynomial[] {
  const n = matrix.length
  const parts = matrix.map((row) => row.map((a) => splitFast(a, fast)))
  const fastEdge = (i: number, j: number): Polynomial => parts[j]![i]!.fast // i → j

  // Pools: states joined by fast reactions
  const pool = new Array<number>(n).fill(-1)
  const pools: number[][] = []
  for (let s = 0; s < n; s++) {
    if (pool[s]! >= 0) continue
    const members = [s]
    pool[s] = pools.length
    for (let k = 0; k < members.length; k++) {
      const i = members[k]!
      for (let j = 0; j < n; j++)
        if (j !== i && pool[j]! < 0 && (!fastEdge(i, j).isZero() || !fastEdge(j, i).isZero())) {
          pool[j] = pools.length
          members.push(j)
        }
    }
    pools.push(members)
  }

  // Within a pool, equilibrium proportions along a spanning tree from its first state
  const proportion = new Array<Fraction | undefined>(n)
  for (const members of pools) {
    proportion[members[0]!] = F.fromPolynomial(Polynomial.one)
    const queue = [members[0]!]
    while (queue.length > 0) {
      const i = queue.shift()!
      for (const j of members) {
        const forward = fastEdge(i, j)
        const reverse = fastEdge(j, i)
        if (forward.isZero() && reverse.isZero()) continue
        if (forward.isZero() || reverse.isZero())
          fail('not-thermodynamically-consistent', 'A fast reaction that only runs one way has no equilibrium')
        const expected = F.mul(proportion[i]!, F.fraction(forward, reverse))
        const known = proportion[j]
        if (known === undefined) {
          proportion[j] = expected
          queue.push(j)
        } else if (!F.sub(known, expected).num.isZero())
          fail('not-thermodynamically-consistent', 'The fast reactions’ equilibria disagree around a loop')
      }
    }
  }

  // Clear each pool's denominators: n_i polynomials in the same proportions
  const numerators = new Array<Polynomial>(n)
  for (const members of pools) {
    let common = Polynomial.one
    for (const i of members) {
      const d = proportion[i]!.den
      if (common.divideExact(d) === undefined) common = common.mul(d)
    }
    let ns = members.map((i) => proportion[i]!.num.mul(common.divideExact(proportion[i]!.den)!))
    let m: Monomial | undefined
    for (const p of ns) m = m === undefined ? p.content().monomial : gcdMonomials(m, p.content().monomial)
    ns = ns.map((p) => p.divideTerm(m ?? []))
    members.forEach((i, k) => (numerators[i] = ns[k]!))
  }

  // The pooled cycle: B[P][R] = Σ_{j∈P, i∈R} slow(A[j][i]) · n_i, whose unknowns are pool amounts / Σ n
  const B = pools.map((P) => pools.map((R) => P.reduce((sum, j) => R.reduce((s, i) => s.add(parts[j]![i]!.slow.mul(numerators[i]!)), sum), Polynomial.zero)))
  const last = pools.length - 1
  const poolWeights = pools.map((_, R) => (pools.length === 1 ? Polynomial.one : cofactor(B, last, R)))
  return Array.from({ length: n }, (_, i) => numerators[i]!.mul(poolWeights[pool[i]!]!))
}
