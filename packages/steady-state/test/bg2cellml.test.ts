/**
 * End to end on CellML that bg2cellml generates: BG-RDF models converted,
 * their maths read back, and the cycle found, solved, written in the common
 * form, checked and reduced.
 */

import { readdirSync, readFileSync } from 'node:fs'
import { fileURLToPath } from 'node:url'

import { Converter, type Runtime } from '@physiomelinks/bg2cellml'
import { parseMathML, toEquations, toMathML, walk, type Equation } from '@physiomelinks/math-ast'
import { DOMParser } from '@xmldom/xmldom'
import createLibCellML from 'libcellml.js'
import * as oxigraph from 'oxigraph'
import { beforeAll, describe, expect, it } from 'vitest'

import { deriveSteadyState, evaluateAst, reduceModel, suggestSpecs, verify } from '../src/index.js'

const ROOT = fileURLToPath(new URL('../../../', import.meta.url))
const BG_RDF = `${ROOT}bg-rdf/`
const FIXTURES = fileURLToPath(new URL('./fixtures/', import.meta.url))

let converter: Converter

beforeAll(async () => {
  const libcellml = await createLibCellML()
  converter = new Converter({ libcellml, oxigraph, domParser: new DOMParser() } as unknown as Runtime)
  const templates: Record<string, string> = {}
  for (const file of readdirSync(`${BG_RDF}templates`)) if (file.endsWith('.ttl')) templates[file] = readFileSync(`${BG_RDF}templates/${file}`, 'utf8')
  expect(converter.loadFramework({ ontology: readFileSync(`${BG_RDF}schema/ontology.ttl`, 'utf8'), templates })).toEqual([])
})

/** Convert a fixture, returning its equations and its variables' initial values. */
async function convert(name: string): Promise<{ equations: Equation[]; values: Record<string, number>; cellml: string }> {
  const result = await converter.convert(readFileSync(`${FIXTURES}${name}`, 'utf8'), { sourceUri: `https://example.org/${name}`, rdfSource: true })
  expect(result.errors).toEqual([])
  const cellml = result.cellml!
  const equations = (cellml.match(/<math[\s\S]*?<\/math>/g) ?? []).flatMap((math) => toEquations(parseMathML(math, { domParser: new DOMParser() as never })))
  const values: Record<string, number> = {}
  for (const m of cellml.matchAll(/<variable\b[^>]*>/g)) {
    const name_ = /\bname="([^"]+)"/.exec(m[0])?.[1]
    const initial = /\binitial_value="([^"]+)"/.exec(m[0])?.[1]
    if (name_ !== undefined && initial !== undefined && !Number.isNaN(Number(initial))) values[name_] = Number(initial)
  }
  return { equations, values, cellml }
}

describe('on bg2cellml’s CellML', () => {
  it('finds, solves and reduces GLUT2, giving Table 3’s km values', async () => {
    const { equations, values } = await convert('glut2.ttl')
    const [cycle, ...others] = suggestSpecs(equations)
    expect(others).toEqual([])
    expect(cycle!.spec.cycleStates).toHaveLength(4)
    expect(cycle!.spec.fast).toHaveLength(2)
    const spec = { ...cycle!.spec, label: 'GLUT2' }
    const result = deriveSteadyState(equations, spec)
    expect(result.issues).toEqual([])
    const common = result.common!
    expect(common.bracket).toBe('difference')
    // km¹, km², km³ (eq 20) from the fixture's Table 3 values
    const lumped = new Map([...common.variables, ...common.definitions].map(({ name, value }) => [name, value]))
    const ks = common.phi.terms.flatMap((t) => (t.k === undefined ? [] : [evaluateAst(lumped.get(t.k.name)!, { ...values, q_tot: 1 }, lumped)])).sort((a, b) => a - b)
    expect(ks[0]).toBeCloseTo(1.4667, 3)
    expect(ks[1]).toBeCloseTo(21.603, 2)
    expect(ks[2]).toBeCloseTo(235.08, 1)
    expect(verify(equations, spec, result, { fastScale: 1e9 }).maxRelativeError).toBeLessThan(1e-5)

    // The reduced model writes back to valid MathML, with no atoms left in it
    const reduced = reduceModel(equations, spec, result)
    const math = toMathML(reduced.equations, { cnUnits: () => 'dimensionless' })
    expect(math).not.toContain('$exp')
    expect(reduced.parameters).toEqual(['q_tot'])
    for (const q of spec.cycleStates) expect(reduced.equations.some((e) => e.left.type === 'Derivative' && [...walk(e.left)].some((n) => n.type === 'Identifier' && n.name === q))).toBe(false)
  })

  it('finds the charge of an electrogenic uniporter from its transformers', async () => {
    const { equations, values } = await convert('uniporter.ttl')
    const ratios = Object.keys(values).filter((v) => v.startsWith('TF_'))
    const [cycle] = suggestSpecs(equations, { ratios, values: { ...values, F: 96485.33212 } })
    expect(cycle!.spec.relations).toBeDefined()
    expect(cycle!.spec.opaque?.some((v) => v.startsWith('u_'))).toBe(true)
    const result = deriveSteadyState(equations, cycle!.spec)
    expect(result.issues).toEqual([])
    const z = result.common!.exponent?.z
    expect(z === undefined ? undefined : Math.abs(Number(z.n) / Number(z.d))).toBe(1)
    expect(result.common!.reversalPotential).toBeDefined()
    // Physical scales, so exp(F u/RT) stays finite
    const physical = { F: 96485.33212, RT: 2578.73, q_C_m: values.q_C_m!, C_C_m: values.C_C_m! }
    expect(verify(equations, cycle!.spec, result, { fastScale: 1e9, values: physical }).maxRelativeError).toBeLessThan(1e-5)
  })
})
