/**
 * The common form of a transporter's steady-state flux (Hunter et al., "The
 * biophysical basis of intracellular sodium homeostasis", eq 24):
 *
 *     v = κ/Φ · { R · (Q_drive/Q_drive^eq) − Q_eq · P · exp(z F u/RT) }
 *
 * R and P are the held species on each side (as `q̄ = K q`, or as
 * concentrations), the driver ratio is a driving reaction's species such as
 * ATP hydrolysis' (1 when there's none), Q_eq gathers the thermodynamic
 * constants, the exponential carries the membrane potential, and Φ, a sum of
 * terms `m/k_m` in the held species, sets how fast the flux is away from
 * equilibrium. Bracket and Φ are read off the solved flux N/D: N must be one
 * forward and one reverse term (a single thermodynamic cycle), and D's term
 * free of held species (when there is one) normalises Φ to `1 + Σ m/k_m`.
 */

import { ci, cn, divide, eq, fn, minus, plus, times, type AstNode, type Equation } from '@physiomelinks/math-ast'

import type { CycleSolution } from './cycle.js'
import { fail } from './issues.js'
import { divideLaurent, inverse, laurent, laurentAst, multiply, pick, rationalAst, type Laurent } from './laurent.js'
import { fractionToAst, type ExpAtom } from './normalise.js'
import * as F from './poly/fraction.js'
import type { Fraction } from './poly/fraction.js'
import { Polynomial, type Term } from './poly/polynomial.js'
import * as Q from './poly/rational.js'
import type { Rational } from './poly/rational.js'

export interface CommonFormOptions {
  /** A name for the transporter, appended to the lumped parameters' names: `GLUT2` gives `kappa_m_GLUT2`. */
  readonly label?: string
  /** The held species' amounts, such as `q_Ao` and `q_Ai`. */
  readonly solutes: readonly string[]
  /** Each held species' thermodynamic constant, by amount (`q_Ao → K_Ao`), to write `K q` as one `qbar_Ao`. */
  readonly pairs?: ReadonlyMap<string, string>
  /** The held species of a driving reaction such as ATP hydrolysis, kept apart as its own ratio. */
  readonly driver?: readonly string[]
  /**
   * `difference`: `[A]_i − Q_eq [A]_o`, with Φ normalised to `1 + Σ …` (GLUT2,
   * Kir); `ratio`: both terms divided by the reverse side (SGLT1, NKA). By
   * default a difference when each side has at most one species, once.
   */
  readonly bracket?: 'difference' | 'ratio'
  /**
   * Write held species as concentrations: each amount's compartment volume
   * (`q_Ao → V_o`), so `qbar = K V [c]` and K V go into κ, Q_eq and the k's.
   */
  readonly concentrations?: ReadonlyMap<string, string>
  /** The membrane potential, when the exp atoms are in one; found from them if not given. */
  readonly potential?: string
}

export interface Lumped {
  readonly name: string
  readonly value: AstNode
}

export interface PhiTerm {
  /** The held species (and exp atoms) of the term. */
  readonly monomial: Laurent
  /** The term is `monomial / k`, when Φ is normalised... */
  readonly k?: { readonly name: string; readonly value: Fraction }
  /** ...or `coefficient · monomial` when it isn't. */
  readonly coefficient?: Fraction
}

export interface CommonForm {
  readonly bracket: 'difference' | 'ratio'
  readonly concentrations: boolean
  readonly kappa: { readonly name: string; readonly value: Fraction }
  /** The held species on the forward and reverse side, not counting the driver's. */
  readonly reactants: Laurent
  readonly products: Laurent
  /** `undefined` when it is exactly 1. */
  readonly qeq?: { readonly name: string; readonly value: Fraction }
  readonly driver?: {
    readonly reactants: Laurent
    readonly products: Laurent
    /** In concentrations, the driving reaction's own equilibrium constant. */
    readonly qeq?: { readonly name: string; readonly value: Fraction }
  }
  /** The exp atoms on the reverse term, with the net charge `z` when they are one exponent. */
  readonly exponent?: { readonly atoms: Laurent; readonly z?: Rational; readonly potential?: string }
  readonly phi: { readonly name: string; readonly prefactor: Laurent; readonly terms: readonly PhiTerm[] }
  /** The membrane potential at which the flux is zero, when there's one exponent. */
  readonly reversalPotential?: AstNode
  /** The grouped or concentration variables used: `qbar_Ao = K_Ao q_Ao`. */
  readonly variables: readonly Lumped[]
  /** κ, Q_eq, the k's and Φ, each as an expression in the model's own parameters. */
  readonly definitions: readonly Lumped[]
  /** `v = κ/Φ · {…}`, in the named lumped parameters. */
  readonly flux: AstNode
  /** The exp atoms the expressions use. */
  readonly atoms: ReadonlyMap<string, ExpAtom>
}

const short = (q: string): string => q.replace(/^(qbar_|q_|c_)/, '')

function renameProducts(p: Polynomial, rename: (m: Term['monomial']) => Term['monomial']): Polynomial {
  return Polynomial.fromTerms(p.terms().map((t) => ({ monomial: laurent(rename(t.monomial)), coefficient: t.coefficient })))
}

/** Rewrite a steady-state flux in the common form, or throw a `SteadyStateError` saying why it can't be. */
export function commonForm(solution: CycleSolution, options: CommonFormOptions): CommonForm {
  const { atoms } = solution
  const suffix = options.label ? `_${options.label}` : ''
  let N = solution.flux.num
  let D = solution.flux.den
  const variables: Lumped[] = []

  // Group K q as one q̄ wherever every term has them to the same power
  const held = new Map<string, string>() // name in the polynomials → model amount
  for (const q of options.solutes) {
    const K = options.pairs?.get(q)
    const pairs = (p: Polynomial) => p.terms().every((t) => (t.monomial.find(([v]) => v === q)?.[1] ?? 0) === (t.monomial.find(([v]) => v === K)?.[1] ?? 0))
    if (K !== undefined && pairs(N) && pairs(D)) {
      const name = `qbar_${short(q)}`
      const rename = (m: Term['monomial']) => m.flatMap(([v, e]): [string, number][] => (v === K ? [] : v === q ? [[name, e]] : [[v, e]]))
      N = renameProducts(N, rename)
      D = renameProducts(D, rename)
      held.set(name, q)
      variables.push({ name, value: times(ci(K), ci(q)) })
    } else held.set(q, q)
  }
  if (held.size === 0)
    fail('no-held-species', 'No held species were found: name the species the cycle moves (constants or states outside it) in `solutes`')
  const isHeld = (v: string) => held.has(v)
  const isAtom = (v: string) => atoms.has(v)
  const drive = [...held.keys(), ...atoms.keys()]

  // The bracket: one forward and one reverse term
  const groups = [...N.coefficientsIn(drive).values()]
  if (groups.length !== 2)
    fail(
      'not-thermodynamically-consistent',
      `The flux's numerator has ${groups.length} terms in the held species, not one forward and one reverse term: is there more than one cycle, such as slippage?`,
    )
  const positive = (p: Polynomial) => Q.sign(p.leadingTerm()!.coefficient) > 0
  const [forward, reverse] = positive(groups[0]!.coefficient) ? [groups[0]!, groups[1]!] : [groups[1]!, groups[0]!]
  if (positive(reverse.coefficient)) fail('not-thermodynamically-consistent', 'Both terms of the flux run the same way, so it has no equilibrium')
  const cf = forward.coefficient
  const cr = reverse.coefficient.neg()
  let qeq = F.fraction(cr, cf)

  const driver = new Set((options.driver ?? []).map((q) => [...held].find(([, model]) => model === q)?.[0] ?? q))
  const mf: Laurent = forward.monomial
  const mr: Laurent = reverse.monomial
  const R = pick(mf, isHeld)
  const P = pick(mr, isHeld)
  const Af = pick(mf, isAtom)
  const X = divideLaurent(pick(mr, isAtom), Af)
  const Rd = pick(R, (v) => driver.has(v))
  const Pd = pick(P, (v) => driver.has(v))
  let R1 = pick(R, (v) => !driver.has(v))
  let P1 = pick(P, (v) => !driver.has(v))
  const simple = (m: Laurent) => m.length <= 1 && m.every(([, e]) => e === 1)
  const form = options.bracket ?? (simple(R1) && simple(P1) ? 'difference' : 'ratio')
  const M = multiply(multiply(Af, Pd), form === 'ratio' ? P1 : [])

  // Φ = (1/M) Σ (D_m/D_0) m, normalised by D's term free of held species when there is one
  // Φ reads 1 + q̄o/k + q̄i/k + q̄o q̄i/k: the constant first, then by degree
  const degree = (m: Laurent) => m.filter(([v]) => !isAtom(v)).reduce((sum, [, e]) => sum + Math.abs(e), 0)
  const dGroups = [...D.coefficientsIn(drive).values()].sort((a, b) => degree(a.monomial) - degree(b.monomial))
  const D0 = dGroups.find((g) => g.monomial.length === 0)?.coefficient
  let kappa = F.fraction(cf, D0 ?? Polynomial.one)
  let phiTerms: { monomial: Laurent; k?: Fraction; coefficient?: Fraction }[] = dGroups.map(({ monomial, coefficient }) =>
    D0 === undefined
      ? { monomial, coefficient: F.fromPolynomial(coefficient) }
      : monomial.length === 0
        ? { monomial }
        : { monomial, k: F.fraction(D0, coefficient) },
  )
  let prefactor = inverse(M)

  // In concentrations: q̄ = (K V)·c, with the factors moved into κ, Q_eq and the k's
  let driverEq: Fraction | undefined
  let Rd1 = Rd
  let Pd1 = Pd
  if (options.concentrations !== undefined) {
    const factorOf = new Map<string, { factor: Polynomial; name: string }>()
    for (const [name, q] of held) {
      const V = options.concentrations.get(q)
      if (V === undefined) continue
      const K = name.startsWith('qbar_') ? options.pairs?.get(q) : undefined
      const factor = K === undefined ? Polynomial.variable(V) : Polynomial.variable(K).mul(Polynomial.variable(V))
      const c = `c_${short(q)}`
      factorOf.set(name, { factor, name: c })
      variables.push({ name: c, value: divide(ci(q), ci(V)) })
    }
    const fac = (m: Laurent): Fraction =>
      m.reduce((f, [v, e]) => {
        const x = factorOf.get(v)
        return x === undefined ? f : F.mul(f, F.pow(F.fromPolynomial(x.factor), e))
      }, F.fromPolynomial(Polynomial.one))
    const rename = (m: Laurent): Laurent => laurent(m.map(([v, e]) => [factorOf.get(v)?.name ?? v, e] as const))
    kappa = F.mul(kappa, form === 'ratio' ? F.div(fac(R1), fac(P1)) : fac(R1))
    qeq = F.mul(qeq, F.div(fac(P1), fac(R1)))
    if (Rd.length + Pd.length > 0) driverEq = F.div(fac(Pd), fac(Rd))
    const facM = fac(M)
    phiTerms = phiTerms.map((t) => {
      const ratio = F.div(fac(t.monomial), facM)
      const monomial = rename(t.monomial)
      if (t.k !== undefined) return { monomial, k: F.div(t.k, ratio) }
      if (t.coefficient !== undefined) return { monomial, coefficient: F.mul(t.coefficient, ratio) }
      // the normalising term: its coefficient is 1/fac(M), which is 1 unless M has held species
      const one = ratio
      return F.isPolynomial(one) && one.num.constantValue() !== undefined && Q.isOne(one.num.constantValue()!)
        ? { monomial }
        : { monomial, coefficient: one }
    })
    prefactor = rename(prefactor)
    R1 = rename(R1)
    P1 = rename(P1)
    Rd1 = rename(Rd)
    Pd1 = rename(Pd)
  }

  // Names
  const names = new Set<string>()
  const unique = (base: string): string => {
    let name = base
    for (let i = 2; names.has(name); i++) name = `${base}_${i}`
    names.add(name)
    return name
  }
  const kappaName = unique(`kappa_m${suffix}`)
  const qeqIsOne = F.isPolynomial(qeq) && qeq.num.constantValue() !== undefined && Q.isOne(qeq.num.constantValue()!)
  const qeqName = qeqIsOne ? undefined : unique(`Q_eq${suffix}`)
  const driverEqName = driverEq === undefined ? undefined : unique(`Q_eq_drive${suffix}`)
  const phiName = unique(`Phi${suffix}`)
  const namedTerms: PhiTerm[] = phiTerms.map((t) => {
    if (t.k === undefined) return t.coefficient === undefined ? { monomial: t.monomial } : { monomial: t.monomial, coefficient: t.coefficient }
    const species = t.monomial.filter(([v]) => !isAtom(v)).map(([v]) => short(v))
    return { monomial: t.monomial, k: { name: unique(`k_${species.join('_') || '0'}${suffix}`), value: t.k } }
  })

  // The exponent, its charge and the reversal potential
  let exponent: CommonForm['exponent']
  let reversalPotential: AstNode | undefined
  if (X.length > 0) {
    const only = X.length === 1 ? atoms.get(X[0]![0]) : undefined
    const z = only === undefined ? undefined : Q.rational(X[0]![1], only.denominator)
    const candidates = only?.base.filter(([v, e]) => e === 1 && !['F', 'RT', 'R', 'T'].includes(v)).map(([v]) => v) ?? []
    const potential = options.potential ?? (candidates.length === 1 ? candidates[0] : undefined)
    exponent = { atoms: X, ...(z === undefined ? {} : { z }), ...(potential === undefined ? {} : { potential }) }
    if (only !== undefined && z !== undefined && potential !== undefined && only.base.some(([v, e]) => v === potential && e === 1)) {
      // R·(Rd/Pd)/(Q_eq·P) = exp(z·base), base = (base/u)·u
      const driverRatio = Rd1.length + Pd1.length > 0 ? divide(laurentAst(Rd1), laurentAst(Pd1)) : undefined
      const top = driverRatio === undefined ? laurentAst(R1) : times(laurentAst(R1), driverEqName === undefined ? driverRatio : divide(driverRatio, ci(driverEqName)))
      const bottom = qeqName === undefined ? laurentAst(P1) : times(ci(qeqName), laurentAst(P1))
      const scale = divideLaurent(only.base, [[potential, 1]])
      reversalPotential = divide(fn('ln', divide(top, bottom)), times(rationalAst(z), laurentAst(scale)))
    }
  }

  // The flux and Φ as expressions
  const driverAst = Rd1.length + Pd1.length > 0 ? (driverEqName === undefined ? divide(laurentAst(Rd1), laurentAst(Pd1)) : divide(divide(laurentAst(Rd1), laurentAst(Pd1)), ci(driverEqName))) : undefined
  const sideAst = (m: Laurent) => (m.length === 0 ? undefined : laurentAst(m))
  const forwardAst = (() => {
    const base = form === 'ratio' ? (R1.length + P1.length === 0 ? undefined : laurentAst(divideLaurent(R1, P1))) : sideAst(R1)
    const parts = [base, driverAst].filter((x): x is AstNode => x !== undefined)
    return parts.length === 0 ? cn(1) : parts.length === 1 ? parts[0]! : times(...parts)
  })()
  const reverseAst = (() => {
    const parts = [qeqName === undefined ? undefined : ci(qeqName), form === 'ratio' ? undefined : sideAst(P1), X.length === 0 ? undefined : laurentAst(X, atoms)].filter(
      (x): x is AstNode => x !== undefined,
    )
    return parts.length === 0 ? cn(1) : parts.length === 1 ? parts[0]! : times(...parts)
  })()
  const flux = times(divide(ci(kappaName), ci(phiName)), minus(forwardAst, reverseAst))

  const termAst = (t: PhiTerm): AstNode => {
    const m = t.monomial.length === 0 ? cn(1) : laurentAst(t.monomial, atoms)
    if (t.k !== undefined) return t.monomial.length === 0 ? divide(cn(1), ci(t.k.name)) : divide(m, ci(t.k.name))
    if (t.coefficient !== undefined) return t.monomial.length === 0 ? fractionToAst(t.coefficient) : times(fractionToAst(t.coefficient), m)
    return m
  }
  const sum = namedTerms.length === 1 ? termAst(namedTerms[0]!) : plus(...namedTerms.map(termAst))
  const phiAst = prefactor.length === 0 ? sum : times(laurentAst(prefactor, atoms), sum)

  const definitions: Lumped[] = [
    { name: kappaName, value: fractionToAst(kappa) },
    ...(qeqName === undefined ? [] : [{ name: qeqName, value: fractionToAst(qeq) }]),
    ...(driverEqName === undefined ? [] : [{ name: driverEqName, value: fractionToAst(driverEq!) }]),
    ...namedTerms.flatMap((t) => (t.k === undefined ? [] : [{ name: t.k.name, value: fractionToAst(t.k.value) }])),
    { name: phiName, value: phiAst },
  ]

  return {
    bracket: form,
    concentrations: options.concentrations !== undefined,
    kappa: { name: kappaName, value: kappa },
    reactants: R1,
    products: P1,
    ...(qeqName === undefined ? {} : { qeq: { name: qeqName, value: qeq } }),
    ...(Rd.length + Pd.length > 0
      ? { driver: { reactants: Rd1, products: Pd1, ...(driverEqName === undefined ? {} : { qeq: { name: driverEqName, value: driverEq! } }) } }
      : {}),
    ...(exponent === undefined ? {} : { exponent }),
    phi: { name: phiName, prefactor, terms: namedTerms },
    ...(reversalPotential === undefined ? {} : { reversalPotential }),
    variables,
    definitions,
    flux,
    atoms,
  }
}

/** The common form's equations: the variables, the lumped parameters, then `v = …`. */
export function commonFormEquations(form: CommonForm, fluxName: string): Equation[] {
  return [...form.variables, ...form.definitions].map(({ name, value }) => eq(ci(name), value)).concat([eq(ci(fluxName), form.flux)])
}
