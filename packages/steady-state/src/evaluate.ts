/**
 * Numeric evaluation of an expression tree, inlining definitions: for
 * checking derived expressions against the model they came from.
 */

import type { AstNode } from '@physiomelinks/math-ast'

export type Definitions = ReadonlyMap<string, AstNode>

const FUNCTIONS: Record<string, (...args: number[]) => number> = {
  exp: Math.exp,
  ln: Math.log,
  log: (x, base = 10) => Math.log(x) / Math.log(base),
  sqrt: Math.sqrt,
  sin: Math.sin,
  cos: Math.cos,
  tan: Math.tan,
  sinh: Math.sinh,
  cosh: Math.cosh,
  tanh: Math.tanh,
  min: Math.min,
  max: Math.max,
  floor: Math.floor,
  ceiling: Math.ceil,
}

/** The value of `node`, with `values` for variables and `definitions` for computed ones. */
export function evaluateAst(node: AstNode, values: ReadonlyMap<string, number> | Readonly<Record<string, number>>, definitions: Definitions = new Map()): number {
  const lookup = values instanceof Map ? values : new Map(Object.entries(values))
  const cache = new Map<string, number>()
  const visiting = new Set<string>()
  const ev = (n: AstNode): number => {
    switch (n.type) {
      case 'Number':
        return n.value
      case 'Identifier': {
        const given = lookup.get(n.name)
        if (given !== undefined) return given
        const cached = cache.get(n.name)
        if (cached !== undefined) return cached
        const definition = definitions.get(n.name)
        if (definition === undefined) throw new Error(`No value for ${n.name}`)
        if (visiting.has(n.name)) throw new Error(`${n.name} is defined in terms of itself`)
        visiting.add(n.name)
        const value = ev(definition)
        visiting.delete(n.name)
        cache.set(n.name, value)
        return value
      }
      case 'Constant':
        if (n.name === 'pi') return Math.PI
        if (n.name === 'exponentiale') return Math.E
        if (n.name === 'infinity') return Infinity
        return NaN
      case 'Add':
        return n.children.reduce((s, c) => s + ev(c), 0)
      case 'Multiply':
        return n.children.reduce((s, c) => s * ev(c), 1)
      case 'Subtract':
        return ev(n.minuend) - ev(n.subtrahend)
      case 'Negate':
        return -ev(n.value)
      case 'Divide':
        return ev(n.numerator) / ev(n.denominator)
      case 'Power':
        return ev(n.base) ** ev(n.exponent)
      case 'Abs':
        return Math.abs(ev(n.value))
      case 'Root':
        return ev(n.radicand) ** (1 / (n.degree === null ? 2 : ev(n.degree)))
      case 'FunctionCall': {
        const f = FUNCTIONS[n.name]
        if (f === undefined) throw new Error(`Can't evaluate ${n.name}`)
        return f(...n.args.map(ev))
      }
      default:
        throw new Error(`Can't evaluate a ${n.type}`)
    }
  }
  return ev(node)
}
