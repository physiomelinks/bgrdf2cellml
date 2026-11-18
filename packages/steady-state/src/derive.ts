/**
 * The package's entry point: a model's equations and what to reduce, to the
 * steady-state flux in the common transporter form.
 */

import { ci, eq, walk, type AstNode, type Equation } from '@physiomelinks/math-ast'

import { solveCycle, splitEquations, type CycleSolution, type CycleSpec } from './cycle.js'
import { SteadyStateError, type Issue } from './issues.js'
import { commonForm, type CommonForm, type CommonFormOptions } from './normalForm.js'
import { fractionToAst } from './normalise.js'

export interface SteadyStateSpec extends CycleSpec {
  /** A name for the transporter, for the lumped parameters' names. */
  readonly label?: string
  /** The held species' amounts; by default every state outside the cycle the flux depends on. */
  readonly solutes?: readonly string[]
  /** Held species of a driving reaction, such as ATP, ADP, Pi and H⁺. */
  readonly driver?: readonly string[]
  readonly bracket?: CommonFormOptions['bracket']
  /** Each held species' compartment volume, to write the result in concentrations. */
  readonly concentrations?: Readonly<Record<string, string>>
  /** The membrane potential, if it can't be found from the exponents. */
  readonly potential?: string
  /**
   * Extra relations between parameters, applied before solving: Kir's
   * charge balance `z4 = −1 + z1 − z2 + z3 + z5 − z6`, say.
   */
  readonly relations?: Readonly<Record<string, AstNode>>
}

export interface SteadyStateResult {
  /** The cycle could be solved. */
  readonly ok: boolean
  /** Why it couldn't be solved, or why the common form couldn't be made (with `ok` still true). */
  readonly issues: readonly Issue[]
  readonly solution?: CycleSolution
  /** The flux as one expression in the model's own variables. */
  readonly flux?: AstNode
  readonly common?: CommonForm
  /** The held species used. */
  readonly solutes: readonly string[]
}

/** Each held species' thermodynamic constant, from chemical potentials written `RT ln(K q)`. */
export function boltzmannPairs(definitions: ReadonlyMap<string, AstNode>, solutes: readonly string[]): Map<string, string> {
  const wanted = new Set(solutes)
  const pairs = new Map<string, string>()
  for (const definition of definitions.values())
    for (const node of walk(definition)) {
      if (node.type !== 'FunctionCall' || node.name !== 'ln') continue
      const arg = node.args[0]
      if (arg?.type !== 'Multiply' || arg.children.length !== 2) continue
      const [a, b] = arg.children
      if (a?.type !== 'Identifier' || b?.type !== 'Identifier') continue
      if (wanted.has(b.name)) pairs.set(b.name, a.name)
      else if (wanted.has(a.name)) pairs.set(a.name, b.name)
    }
  return pairs
}

export function deriveSteadyState(equations: readonly Equation[], spec: SteadyStateSpec): SteadyStateResult {
  const withRelations = [...equations, ...Object.entries(spec.relations ?? {}).map(([name, value]) => eq(ci(name), value))]
  const solved = solveCycle(withRelations, spec)
  if (!solved.ok) return { ok: false, issues: solved.issues, solutes: [] }
  const solution = solved.value
  const { definitions, rates } = splitEquations(withRelations)
  const cycle = new Set(spec.cycleStates)
  const used = new Set([...solution.flux.num.variables(), ...solution.flux.den.variables(), ...[...solution.atoms.values()].flatMap((a) => a.base.map(([v]) => v))])
  const solutes = spec.solutes ?? [...rates.keys()].filter((q) => !cycle.has(q) && used.has(q))
  const flux = fractionToAst(solution.flux, solution.atoms)
  try {
    const common = commonForm(solution, {
      solutes,
      pairs: boltzmannPairs(definitions, solutes),
      ...(spec.label === undefined ? {} : { label: spec.label }),
      ...(spec.driver === undefined ? {} : { driver: spec.driver }),
      ...(spec.bracket === undefined ? {} : { bracket: spec.bracket }),
      ...(spec.concentrations === undefined ? {} : { concentrations: new Map(Object.entries(spec.concentrations)) }),
      ...(spec.potential === undefined ? {} : { potential: spec.potential }),
    })
    return { ok: true, issues: [], solution, flux, common, solutes }
  } catch (e) {
    if (e instanceof SteadyStateError) return { ok: true, issues: [e.issue], solution, flux, solutes }
    throw e
  }
}
