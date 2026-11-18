/**
 * Why a model can't be reduced, in the shape bg2cellml's issues have: a code
 * a host can switch on and a message a person can read.
 */

export type IssueCode =
  | 'unsupported-construct'
  | 'unsupported-logarithm'
  | 'nonlinear-in-states'
  | 'no-conservation'
  | 'disconnected-cycle'
  | 'unknown-variable'
  | 'circular-definition'
  | 'not-thermodynamically-consistent'

export interface Issue {
  readonly code: IssueCode
  readonly message: string
  /** The variable the issue is about, when there is one. */
  readonly variable?: string
}

/** Thrown inside the algebra and turned into an issue at the API boundary. */
export class SteadyStateError extends Error {
  readonly issue: Issue
  constructor(issue: Issue) {
    super(issue.message)
    this.name = 'SteadyStateError'
    this.issue = issue
  }
}

export const fail = (code: IssueCode, message: string, variable?: string): never => {
  throw new SteadyStateError(variable === undefined ? { code, message } : { code, message, variable })
}
