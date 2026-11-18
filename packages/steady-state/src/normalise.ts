/**
 * From a model's equations to exact quotients of polynomials.
 *
 * Bond-graph CellML writes reaction rates as `κ(exp(u_f/RT) − exp(u_r/RT))`
 * with chemical potentials `u = RT ln(K q)`. Inlined, an exponent becomes a
 * sum of terms, each a rational times a product of variables (`RT ln(K q)/RT`)
 * or a logarithm. Terms whose variable part cancels to 1 with a logarithm give
 * powers, `exp(2 ln(K q)) = (K q)²`; every other term, such as `−z F u/RT`
 * for a membrane potential, becomes an *exp atom*, a new variable standing for
 * `exp(base/L)`, so that `exp(a + b) = exp(a)·exp(b)` and the result is a
 * polynomial in the model's variables and the atoms.
 *
 * Each atom's L is the least common denominator of every coefficient its base
 * gets anywhere in the expressions given, so `exp(F u/(2RT))` and
 * `exp(F u/RT)` are one atom, to the powers 1 and 2. That's why
 * `normalise()` takes all the expressions at once.
 */

import { cn, ci, divide, exp, minus, plus, power, times, type AstNode } from '@physiomelinks/math-ast'

import { fail } from './issues.js'
import * as F from './poly/fraction.js'
import type { Fraction } from './poly/fraction.js'
import { Polynomial } from './poly/polynomial.js'
import * as Q from './poly/rational.js'
import type { Rational } from './poly/rational.js'

/** A product of variables to integer powers, which may be negative: sorted by name, no zero powers. */
export type Laurent = readonly (readonly [string, number])[]

/** A variable standing for `exp(base/denominator)`. */
export interface ExpAtom {
  /** Its name in the polynomials; never a CellML identifier (it starts with `$`). */
  readonly name: string
  readonly base: Laurent
  readonly denominator: number
}

export interface NormaliseOptions {
  /** Each computed variable's defining expression (the right side of `name = …`). */
  readonly definitions: ReadonlyMap<string, AstNode>
  /** Variables, such as a cycle's states, that anything depending on must be inlined. */
  readonly states?: Iterable<string>
  /**
   * Variables never to inline, kept as they are: a membrane potential, say,
   * so an exponent is written in it rather than in the charge it comes from.
   * Outside exponents, a definition is inlined only when it depends on a
   * state or has a logarithm; inside them, always (unless it's opaque).
   */
  readonly opaque?: Iterable<string>
  /** Variables always to inline. */
  readonly inline?: Iterable<string>
}

export interface Normalised {
  /** Each expression given, as a quotient of polynomials in variables and atoms. */
  readonly fractions: readonly Fraction[]
  /** The exp atoms the fractions use, by name. */
  readonly atoms: ReadonlyMap<string, ExpAtom>
}

export const ATOM_PREFIX = '$exp'
export const isAtom = (name: string): boolean => name.startsWith(ATOM_PREFIX)

const laurentKey = (m: Laurent): string => (m.length === 0 ? '1' : m.map(([v, e]) => (e === 1 ? v : `${v}^${e}`)).join('*'))

function multiplyLaurent(a: Laurent, b: Laurent, sign = 1): Laurent {
  const powers = new Map(a)
  for (const [v, e] of b) powers.set(v, (powers.get(v) ?? 0) + sign * e)
  return [...powers].filter(([, e]) => e !== 0).sort(([x], [y]) => (x < y ? -1 : x > y ? 1 : 0))
}

/** A term of an exponent: `c · m`, times `ln(log)` when there is one. */
interface LogTerm {
  readonly c: Rational
  readonly m: Laurent
  readonly log: Fraction | null
}

const lcm = (a: number, b: number): number => {
  const g = (x: number, y: number): number => (y === 0 ? x : g(y, x % y))
  return (a / g(a, b)) * b
}

/** The exact value of a `<cn>`, as written when it was written in e-notation. */
function numberValue(node: Extract<AstNode, { type: 'Number' }>): Rational {
  if (node.scientific !== undefined) return Q.fromDecimal(`${node.scientific.mantissa}e${node.scientific.exponent}`)
  return Q.fromNumber(node.value)
}

/**
 * Turn expressions into quotients of polynomials, inlining definitions as
 * `options` says. Throws a `SteadyStateError` for what it can't represent:
 * piecewise, trigonometric and other functions, a logarithm outside an
 * exponent, or a non-integer power of a logarithm's argument.
 */
export function normalise(expressions: readonly AstNode[], options: NormaliseOptions): Normalised {
  const { definitions } = options
  const states = new Set(options.states ?? [])
  const opaque = new Set(options.opaque ?? [])
  const forced = new Set(options.inline ?? [])

  // Whether a definition depends on a state, or has a logarithm, through other definitions
  const dependsOnState = new Map<string, boolean>()
  const hasLog = new Map<string, boolean>()
  const visiting = new Set<string>()
  const scan = (name: string): void => {
    if (dependsOnState.has(name)) return
    const definition = definitions.get(name)
    if (definition === undefined) {
      dependsOnState.set(name, states.has(name))
      hasLog.set(name, false)
      return
    }
    if (visiting.has(name)) fail('circular-definition', `${name} is defined in terms of itself`, name)
    visiting.add(name)
    let state = false
    let log = false
    const walkNode = (n: AstNode): void => {
      if (n.type === 'Identifier') {
        scan(n.name)
        state ||= states.has(n.name) || dependsOnState.get(n.name)!
        log ||= hasLog.get(n.name)!
        return
      }
      if (n.type === 'FunctionCall' && (n.name === 'ln' || n.name === 'log')) log = true
      for (const child of childNodes(n)) walkNode(child)
    }
    walkNode(definition)
    visiting.delete(name)
    dependsOnState.set(name, state)
    hasLog.set(name, log)
  }
  const inlined = (name: string): AstNode | undefined => {
    const definition = definitions.get(name)
    if (definition === undefined || opaque.has(name)) return undefined
    if (forced.has(name)) return definition
    scan(name)
    return dependsOnState.get(name) || hasLog.get(name) ? definition : undefined
  }

  // Pass 1 records each atom base's coefficient denominators; pass 2 builds the polynomials
  const denominators = new Map<string, { base: Laurent; denominator: number }>()
  let atoms: Map<string, ExpAtom & { key: string }> | undefined

  const atomPower = (base: Laurent, c: Rational): Fraction => {
    const key = laurentKey(base)
    if (atoms === undefined) {
      const entry = denominators.get(key) ?? { base, denominator: 1 }
      entry.denominator = lcm(entry.denominator, Number(c.d))
      denominators.set(key, entry)
      return F.fromPolynomial(Polynomial.one)
    }
    const atom = atoms.get(key)!
    const p = Q.mul(c, Q.rational(atom.denominator))
    const k = Number(p.n) // p.d is 1 by construction of the denominator
    const x = Polynomial.variable(atom.name, Math.abs(k))
    return k >= 0 ? F.fromPolynomial(x) : F.fraction(Polynomial.one, x)
  }

  const logProduct = (a: readonly LogTerm[], b: readonly LogTerm[]): LogTerm[] => {
    const result: LogTerm[] = []
    for (const x of a)
      for (const y of b) {
        if (x.log !== null && y.log !== null) fail('unsupported-logarithm', 'A product of two logarithms in an exponent')
        result.push({ c: Q.mul(x.c, y.c), m: multiplyLaurent(x.m, y.m), log: x.log ?? y.log })
      }
    return collectTerms(result)
  }

  const argument = (node: AstNode): LogTerm[] => {
    switch (node.type) {
      case 'Number':
        return [{ c: numberValue(node), m: [], log: null }]
      case 'Identifier': {
        // In an exponent everything is inlined, so sums split into atoms, unless it's opaque
        const definition = opaque.has(node.name) ? undefined : definitions.get(node.name)
        if (definition !== undefined) {
          if (visiting.has(node.name)) fail('circular-definition', `${node.name} is defined in terms of itself`, node.name)
          visiting.add(node.name)
          const result = argument(definition)
          visiting.delete(node.name)
          return result
        }
        return [{ c: Q.ONE, m: [[node.name, 1]], log: null }]
      }
      case 'Constant':
        if (node.name === 'pi') return [{ c: Q.ONE, m: [['pi', 1]], log: null }]
        return fail('unsupported-construct', `The constant ${node.name} in an exponent`)
      case 'Add':
        return collectTerms(node.children.flatMap(argument))
      case 'Subtract':
        return collectTerms([...argument(node.minuend), ...negateTerms(argument(node.subtrahend))])
      case 'Negate':
        return negateTerms(argument(node.value))
      case 'Multiply':
        return node.children.map(argument).reduce(logProduct, [{ c: Q.ONE, m: [], log: null }])
      case 'Divide': {
        const d = argument(node.denominator)
        const only = d[0]
        if (d.length !== 1 || only === undefined || only.log !== null || Q.isZero(only.c))
          return fail('unsupported-construct', 'An exponent divided by a sum or a logarithm')
        return logProduct(argument(node.numerator), [{ c: Q.div(Q.ONE, only.c), m: multiplyLaurent([], only.m, -1), log: null }])
      }
      case 'Power': {
        const k = integerExponent(node.exponent)
        const base = argument(node.base)
        if (k >= 0) {
          let result: LogTerm[] = [{ c: Q.ONE, m: [], log: null }]
          for (let i = 0; i < k; i++) result = logProduct(result, base)
          return result
        }
        const only = base[0]
        if (base.length !== 1 || only === undefined || only.log !== null)
          return fail('unsupported-construct', 'A negative power of a sum in an exponent')
        return [{ c: Q.pow(only.c, k), m: multiplyLaurent([], only.m, k), log: null }]
      }
      case 'FunctionCall':
        if (node.name === 'ln' && node.args.length === 1) return [{ c: Q.ONE, m: [], log: value(node.args[0]!) }]
        return fail('unsupported-construct', `The function ${node.name} in an exponent`)
      default:
        return fail('unsupported-construct', `A ${node.type} in an exponent`)
    }
  }

  const integerExponent = (node: AstNode): number => {
    const v = value(node)
    const constant = v.num.constantValue()
    const k = constant !== undefined && F.isPolynomial(v) ? Q.asInteger(constant) : undefined
    if (k === undefined) return fail('unsupported-construct', 'A power that isn’t a whole number')
    return k
  }

  const expOf = (terms: readonly LogTerm[]): Fraction => {
    let result = F.fromPolynomial(Polynomial.one)
    for (const term of terms) {
      if (term.log !== null) {
        const k = Q.asInteger(term.c)
        if (term.m.length !== 0 || k === undefined)
          fail('unsupported-logarithm', `exp of ${Q.toString(term.c)}${term.m.length ? `·${laurentKey(term.m)}` : ''}·ln(…): only whole-number multiples of a logarithm can be simplified`)
        result = F.mul(result, F.pow(term.log, k!))
      } else result = F.mul(result, atomPower(term.m, term.c))
    }
    return result
  }

  let memo = new Map<string, Fraction>()
  const value = (node: AstNode): Fraction => {
    switch (node.type) {
      case 'Number':
        return F.fromPolynomial(Polynomial.constant(numberValue(node)))
      case 'Identifier': {
        const definition = inlined(node.name)
        if (definition === undefined) return F.fromPolynomial(Polynomial.variable(node.name))
        const cached = memo.get(node.name)
        if (cached !== undefined) return cached
        const result = value(definition)
        memo.set(node.name, result)
        return result
      }
      case 'Constant':
        if (node.name === 'pi') return F.fromPolynomial(Polynomial.variable('pi'))
        if (node.name === 'exponentiale') return atomPower([], Q.ONE)
        return fail('unsupported-construct', `The constant ${node.name}`)
      case 'Add':
        return node.children.map(value).reduce(F.add)
      case 'Subtract':
        return F.sub(value(node.minuend), value(node.subtrahend))
      case 'Negate':
        return F.neg(value(node.value))
      case 'Multiply':
        return node.children.map(value).reduce(F.mul)
      case 'Divide': {
        const d = value(node.denominator)
        if (d.num.isZero()) return fail('unsupported-construct', 'Division by zero')
        return F.div(value(node.numerator), d)
      }
      case 'Power':
        return F.pow(value(node.base), integerExponent(node.exponent))
      case 'FunctionCall':
        if (node.name === 'exp' && node.args.length === 1) return expOf(argument(node.args[0]!))
        if (node.name === 'ln' || node.name === 'log')
          return fail('unsupported-logarithm', 'A logarithm outside an exponent, such as a chemical potential the flux uses directly')
        return fail('unsupported-construct', `The function ${node.name}`)
      default:
        return fail('unsupported-construct', `A ${node.type}`)
    }
  }

  for (const e of expressions) value(e)
  atoms = new Map(
    [...denominators]
      .sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0))
      .map(([key, { base, denominator }], i) => [key, { key, name: `${ATOM_PREFIX}${i + 1}`, base, denominator }]),
  )
  memo = new Map()
  const fractions = expressions.map(value)
  return {
    fractions,
    atoms: new Map([...atoms.values()].map(({ name, base, denominator }) => [name, { name, base, denominator }])),
  }
}

function negateTerms(terms: readonly LogTerm[]): LogTerm[] {
  return terms.map((t) => ({ ...t, c: Q.neg(t.c) }))
}

function collectTerms(terms: readonly LogTerm[]): LogTerm[] {
  const merged = new Map<string, LogTerm>()
  for (const t of terms) {
    const key = `${laurentKey(t.m)}|${t.log === null ? '' : F.toString(t.log)}`
    const existing = merged.get(key)
    merged.set(key, existing === undefined ? t : { ...existing, c: Q.add(existing.c, t.c) })
  }
  return [...merged.values()].filter((t) => !Q.isZero(t.c))
}

function childNodes(n: AstNode): AstNode[] {
  switch (n.type) {
    case 'Add':
    case 'Multiply':
    case 'And':
    case 'Or':
    case 'Xor':
      return [...n.children]
    case 'Subtract':
      return [n.minuend, n.subtrahend]
    case 'Negate':
    case 'Abs':
    case 'Not':
      return [n.value]
    case 'Divide':
      return [n.numerator, n.denominator]
    case 'Power':
      return [n.base, n.exponent]
    case 'Root':
      return n.degree === null ? [n.radicand] : [n.radicand, n.degree]
    case 'FunctionCall':
      return [...n.args]
    case 'Equal':
    case 'Less':
    case 'Greater':
    case 'LessEqual':
    case 'GreaterEqual':
    case 'NotEqual':
      return [n.left, n.right]
    case 'Piecewise':
      return [...n.pieces.flatMap((p) => [p.value, p.condition]), ...(n.otherwise === null ? [] : [n.otherwise])]
    case 'Derivative':
      return [n.expression]
    default:
      return []
  }
}

//==============================================================================
// Back to expression trees

/** Each atom's value at the given values of its base's variables. */
export function atomValues(atoms: ReadonlyMap<string, ExpAtom>, values: ReadonlyMap<string, number> | Readonly<Record<string, number>>): Map<string, number> {
  const lookup = values instanceof Map ? values : new Map(Object.entries(values))
  const result = new Map(lookup)
  for (const atom of atoms.values()) {
    const base = atom.base.reduce((p, [v, e]) => {
      const x = lookup.get(v)
      if (x === undefined) throw new Error(`No value for ${v}`)
      return p * x ** e
    }, 1)
    result.set(atom.name, Math.exp(base / atom.denominator))
  }
  return result
}

function rationalAst(r: Rational): AstNode {
  return r.d === 1n ? cn(Number(r.n)) : divide(cn(Number(r.n)), cn(Number(r.d)))
}

function productAst(factors: AstNode[]): AstNode {
  return factors.length === 0 ? cn(1) : factors.length === 1 ? factors[0]! : times(...factors)
}

function laurentAst(m: Laurent): AstNode {
  const factor = ([v, e]: readonly [string, number]): AstNode => (Math.abs(e) === 1 ? ci(v) : power(ci(v), cn(Math.abs(e))))
  const top = productAst(m.filter(([, e]) => e > 0).map(factor))
  const bottom = m.filter(([, e]) => e < 0).map(factor)
  return bottom.length === 0 ? top : divide(top, productAst(bottom))
}

/** `exp(r · base)` */
function atomAst(atom: ExpAtom, k: number): AstNode {
  const r = Q.rational(k, atom.denominator)
  const base = laurentAst(atom.base)
  const magnitude = Q.sign(r) < 0 ? Q.neg(r) : r
  const scaled = Q.isOne(magnitude) ? base : times(rationalAst(magnitude), base)
  return exp(Q.sign(r) < 0 ? minus(scaled) : scaled)
}

/** A polynomial as an expression tree, with its atoms written as `exp(…)`. */
export function polynomialToAst(p: Polynomial, atoms: ReadonlyMap<string, ExpAtom> = new Map(), { factor = true }: { factor?: boolean } = {}): AstNode {
  if (p.isZero()) return cn(0)
  // What every term shares comes out front: K_2 κ_1 (K_1 + K_4), not K_1 K_2 κ_1 + K_2 K_4 κ_1
  if (factor && p.size > 1) {
    const { monomial: m, coefficient } = p.content()
    if (m.length > 0 || !Q.isOne(coefficient)) {
      const common = polynomialToAst(Polynomial.term(coefficient, m), atoms, { factor: false })
      const rest = polynomialToAst(p.divideTerm(m, coefficient), atoms, { factor: false })
      return common.type === 'Multiply' ? times(...common.children, rest) : times(common, rest)
    }
  }
  const termAst = (coefficient: Rational, m: readonly (readonly [string, number])[]): AstNode => {
    const factors: AstNode[] = []
    if (!Q.isOne(coefficient)) factors.push(rationalAst(coefficient))
    for (const [v, e] of m) {
      const atom = atoms.get(v)
      if (atom !== undefined) factors.push(atomAst(atom, e))
      else factors.push(e === 1 ? ci(v) : power(ci(v), cn(e)))
    }
    return productAst(factors)
  }
  let result: AstNode | undefined
  for (const { monomial: m, coefficient } of p.terms()) {
    const negative = Q.sign(coefficient) < 0
    const term = termAst(negative ? Q.neg(coefficient) : coefficient, m)
    if (result === undefined) result = negative ? minus(term) : term
    else if (negative) result = minus(result, term)
    else result = result.type === 'Add' ? plus(...result.children, term) : plus(result, term)
  }
  return result!
}

/** A fraction as an expression tree. */
export function fractionToAst(f: Fraction, atoms: ReadonlyMap<string, ExpAtom> = new Map()): AstNode {
  const num = polynomialToAst(f.num, atoms)
  return F.isPolynomial(f) ? num : divide(num, polynomialToAst(f.den, atoms))
}

