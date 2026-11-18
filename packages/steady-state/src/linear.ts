/**
 * A rate of change as a linear combination of reaction fluxes: bond-graph
 * CellML writes `dq/dt = v_q`, `v_q = v_R1 − 2 v_R2`, and so on down to
 * each reaction's own rate `v_R = κ(…)`. Expanding through the sums and
 * constant multiples finds the reactions and their coefficients.
 */

import { divide, minus, times, walk, type AstNode } from '@physiomelinks/math-ast'

export interface LinearTerm {
  /** The reaction flux variable, or `undefined` for an expression that isn't one. */
  readonly flux?: string
  /** Its coefficient, free of the dependent variables. */
  readonly coefficient: AstNode
  /** The expression, when it isn't a flux variable. */
  readonly node?: AstNode
}

/** Which variables (transitively through definitions) depend on any of `on`. */
export function dependents(definitions: ReadonlyMap<string, AstNode>, on: Iterable<string>): (name: string) => boolean {
  const roots = new Set(on)
  const memo = new Map<string, boolean>()
  const depends = (name: string, seen = new Set<string>()): boolean => {
    if (roots.has(name)) return true
    const cached = memo.get(name)
    if (cached !== undefined) return cached
    const definition = definitions.get(name)
    if (definition === undefined || seen.has(name)) return false
    seen.add(name)
    let result = false
    for (const n of walk(definition)) if (n.type === 'Identifier' && depends(n.name, seen)) result = true
    memo.set(name, result)
    return result
  }
  return (name) => depends(name)
}

/** Whether an expression depends on the variables `isDependent` picks. */
export function dependsOn(node: AstNode, isDependent: (name: string) => boolean): boolean {
  for (const n of walk(node)) if (n.type === 'Identifier' && isDependent(n.name)) return true
  return false
}

const scale = (c: AstNode, by: AstNode): AstNode => (c.type === 'Number' && c.value === 1 ? by : times(c, by))
const negate = (c: AstNode): AstNode => (c.type === 'Number' ? { ...c, value: -c.value } : minus(c))

/**
 * Expand `node` into terms `coefficient · flux` through sums, differences,
 * negations and multiples by expressions that don't depend on the chosen
 * variables. Terms that don't depend on them at all are dropped.
 */
export function expandLinear(node: AstNode, definitions: ReadonlyMap<string, AstNode>, isDependent: (name: string) => boolean): LinearTerm[] {
  const one: AstNode = { type: 'Number', value: 1 }
  const expand = (n: AstNode, c: AstNode): LinearTerm[] => {
    if (!dependsOn(n, isDependent)) return []
    switch (n.type) {
      case 'Add':
        return n.children.flatMap((child) => expand(child, c))
      case 'Subtract':
        return [...expand(n.minuend, c), ...expand(n.subtrahend, negate(c))]
      case 'Negate':
        return expand(n.value, negate(c))
      case 'Multiply': {
        const dependent = n.children.filter((child) => dependsOn(child, isDependent))
        if (dependent.length === 1) {
          const others = n.children.filter((child) => child !== dependent[0])
          return expand(dependent[0]!, others.reduce(scale, c))
        }
        return [{ coefficient: c, node: n }]
      }
      case 'Divide':
        if (!dependsOn(n.denominator, isDependent)) return expand(n.numerator, scale(c, divide(one, n.denominator)))
        return [{ coefficient: c, node: n }]
      case 'Identifier': {
        const definition = definitions.get(n.name)
        if (definition !== undefined && isDependent(n.name)) {
          const inner = expand(definition, c)
          // A reaction's own rate isn't a combination of others: it is the flux
          if (inner.every((t) => t.flux !== undefined)) return inner
        }
        return [{ flux: n.name, coefficient: c }]
      }
      default:
        return [{ coefficient: c, node: n }]
    }
  }
  return expand(node, one)
}
