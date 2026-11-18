/**
 * The steady-state version of a whole model: the cycle's states and
 * reactions replaced by the derived flux, and everything else kept.
 */

import { ci, eq, mapChildren, minus, times, walk, type AstNode, type Equation } from '@physiomelinks/math-ast'

import { solveCycle, splitEquations } from './cycle.js'
import type { SteadyStateResult, SteadyStateSpec } from './derive.js'
import { dependents, expandLinear } from './linear.js'
import { commonFormEquations } from './normalForm.js'
import { fractionToAst, isAtom } from './normalise.js'
import * as F from './poly/fraction.js'
import type { Fraction } from './poly/fraction.js'
import * as Q from './poly/rational.js'

export interface ReducedSection {
  readonly title: string
  /** The variables the section's equations define, in order. */
  readonly variables: readonly string[]
}

export interface ReducedModel {
  readonly equations: readonly Equation[]
  readonly sections: readonly ReducedSection[]
  /** The steady-state flux's variable. */
  readonly flux: string
  /** Variables the reduced model has that the original doesn't define or give values to, such as the cycle's total. */
  readonly parameters: readonly string[]
}

export interface ReduceOptions {
  /** The flux's variable name; `v_<label>` (or `v_ss`) if not given. */
  readonly flux?: string
  /**
   * Variables of the cycle to keep, at their steady state: a reaction's flux as its multiple of the flux
   * (`v_R1 = v_ss`, `v_Na = 2 v_ss`), so what refers to it by name, such as an annotation or a port, still can.
   */
  readonly keep?: Iterable<string>
}

const leftName = (e: Equation): string =>
  e.left.type === 'Identifier' ? e.left.name : e.left.type === 'Derivative' && e.left.expression.type === 'Identifier' ? e.left.expression.name : ''

/**
 * Replace the cycle in `equations` by its steady state. Every use of a cycle
 * reaction's flux becomes a multiple of the derived flux (or, when it isn't
 * one, its own steady-state expression), every use of a cycle state its
 * steady-state amount, and the cycle's own equations go.
 */
export function reduceModel(equations: readonly Equation[], spec: SteadyStateSpec, result: SteadyStateResult, options: ReduceOptions = {}): ReducedModel {
  if (result.solution === undefined) throw new Error('Nothing to reduce: the cycle wasn’t solved')
  const withRelations = [...equations, ...Object.entries(spec.relations ?? {}).map(([name, value]) => eq(ci(name), value))]
  const { definitions, rates } = splitEquations(withRelations)
  const cycle = new Set(spec.cycleStates)
  const isDependent = dependents(definitions, cycle)
  const main = result.solution.flux
  const atoms = result.solution.atoms
  const total = result.solution.total

  const taken = new Set([...definitions.keys(), ...rates.keys()])
  let flux = options.flux ?? `v_${spec.label ?? 'ss'}`
  for (let i = 2; taken.has(flux); i++) flux = `${options.flux ?? `v_${spec.label ?? 'ss'}`}_${i}`

  // The variables a ratio to the flux may use: not held species, atoms or the total
  const heldOrNew = new Set([...result.solutes, total, ...(result.common?.variables.map((v) => v.name) ?? [])])
  const isConstantRatio = (f: Fraction): boolean => [...f.num.variables(), ...f.den.variables()].every((v) => !heldOrNew.has(v) && !isAtom(v))

  // What each dependent variable becomes
  const replacements = new Map<string, AstNode>()
  const replacement = (name: string): AstNode => {
    const known = replacements.get(name)
    if (known !== undefined) return known
    const solved = solveCycle(withRelations, { ...spec, flux: name })
    let node: AstNode
    if (!solved.ok) node = ci(name)
    else {
      const ratio = F.div(solved.value.flux, main)
      const constant = F.isPolynomial(ratio) ? ratio.num.constantValue() : undefined
      if (constant !== undefined && Q.isOne(constant)) node = ci(flux)
      else if (constant !== undefined && Q.equal(constant, Q.rational(-1))) node = minus(ci(flux))
      else if (isConstantRatio(ratio)) node = times(fractionToAst(ratio), ci(flux))
      else node = fractionToAst(solved.value.flux, atoms)
    }
    replacements.set(name, node)
    return node
  }
  const keep = new Set(options.keep ?? [])
  /** A kept variable: a reaction's flux as its multiple of the flux, a sum of fluxes rewritten. */
  const atomicOrRewritten = (name: string, right: AstNode): AstNode => {
    const terms = expandLinear(ci(name), definitions, isDependent)
    return terms.length === 1 && terms[0]!.flux === name ? replacement(name) : rewrite(right)
  }

  /** Inline sums of fluxes; replace a reaction's flux or a state by its steady state. */
  const rewrite = (node: AstNode): AstNode => {
    if (node.type === 'Identifier') {
      if (!isDependent(node.name) || keep.has(node.name)) return node
      const definition = definitions.get(node.name)
      const terms = definition === undefined ? [] : expandLinear(ci(node.name), definitions, isDependent)
      const atomic = definition === undefined || (terms.length === 1 && terms[0]!.flux === node.name)
      return atomic ? replacement(node.name) : rewrite(definition)
    }
    return mapChildren(node, rewrite)
  }

  // Which definitions the original model uses, to find its outputs (used by nothing)
  const used = new Set<string>()
  for (const e of withRelations) for (const n of walk(e.right)) if (n.type === 'Identifier') used.add(n.name)

  const flux_: Equation[] = result.common !== undefined ? commonFormEquations(result.common, flux) : [eq(ci(flux), fractionToAst(main, atoms))]
  const held: Equation[] = []
  const outputs: Equation[] = []
  for (const e of withRelations) {
    const name = leftName(e)
    if (e.left.type === 'Derivative') {
      if (cycle.has(name)) continue
      if (dependsOnCycle(e.right, isDependent)) held.push(eq(e.left, rewrite(e.right)))
    } else if (e.left.type === 'Identifier' && isDependent(name) && (keep.has(name) || !used.has(name))) {
      outputs.push(eq(e.left, keep.has(name) ? atomicOrRewritten(name, e.right) : rewrite(e.right)))
    }
  }

  // Keep what the kept equations need, and the rest of the model unchanged
  const kept = new Map<string, Equation>()
  const rewritten = new Set([...held, ...outputs].map(leftName))
  for (const e of withRelations) {
    const name = leftName(e)
    if (rewritten.has(name) && !(e.left.type === 'Identifier' && rates.has(name))) continue
    if (e.left.type === 'Derivative' && (cycle.has(name) || rewritten.has(name))) continue
    if (e.left.type === 'Identifier' && isDependent(name)) continue
    kept.set(`${e.left.type}:${name}`, e)
  }
  // Drop unchanged definitions only the cycle needed
  const needed = new Set<string>()
  const need = (node: AstNode) => {
    for (const n of walk(node))
      if (n.type === 'Identifier' && !needed.has(n.name)) {
        needed.add(n.name)
        const d = definitions.get(n.name)
        if (d !== undefined && !isDependent(n.name)) need(d)
      }
  }
  for (const e of [...flux_, ...held, ...outputs]) need(e.right)
  for (const e of kept.values()) if (e.left.type === 'Derivative' || !used.has(leftName(e))) need(e.right)
  const unchanged = [...kept.values()].filter((e) => e.left.type === 'Derivative' || needed.has(leftName(e)) || !used.has(leftName(e)))

  const defined = new Set([...definitions.keys(), ...rates.keys(), ...flux_.map(leftName)])
  const parameters = [...new Set([...flux_, ...held, ...outputs].flatMap((e) => [...walk(e.right)].flatMap((n) => (n.type === 'Identifier' ? [n.name] : []))))].filter(
    (v) => !defined.has(v) && !usedOriginally(v, withRelations),
  )

  return {
    equations: [...flux_, ...held, ...outputs, ...unchanged],
    sections: [
      { title: 'Steady-state flux', variables: flux_.map(leftName) },
      ...(held.length > 0 ? [{ title: 'Held species', variables: held.map(leftName) }] : []),
      ...(outputs.length > 0 ? [{ title: 'Kept at steady state', variables: outputs.map(leftName) }] : []),
      ...(unchanged.length > 0 ? [{ title: 'Unchanged', variables: unchanged.map(leftName) }] : []),
    ],
    flux,
    parameters,
  }
}

function dependsOnCycle(node: AstNode, isDependent: (name: string) => boolean): boolean {
  for (const n of walk(node)) if (n.type === 'Identifier' && isDependent(n.name)) return true
  return false
}

function usedOriginally(name: string, equations: readonly Equation[]): boolean {
  for (const e of equations) for (const n of [...walk(e.left), ...walk(e.right)]) if (n.type === 'Identifier' && n.name === name) return true
  return false
}
