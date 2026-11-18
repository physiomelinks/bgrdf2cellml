/**
 * Determinants and cofactors of square matrices of polynomials, by Bareiss's
 * fraction-free elimination: every division it makes is exact, so entries stay
 * polynomials and never grow into fractions.
 */

import { Polynomial } from './polynomial.js'

export type PolynomialMatrix = readonly (readonly Polynomial[])[]

function assertSquare(matrix: PolynomialMatrix): number {
  const n = matrix.length
  for (const row of matrix) if (row.length !== n) throw new RangeError('The matrix is not square')
  return n
}

/** The determinant; 1 for the empty matrix. */
export function determinant(matrix: PolynomialMatrix): Polynomial {
  const n = assertSquare(matrix)
  if (n === 0) return Polynomial.one
  if (n === 1) return matrix[0]![0]!
  if (n === 2) return matrix[0]![0]!.mul(matrix[1]![1]!).sub(matrix[0]![1]!.mul(matrix[1]![0]!))

  const m = matrix.map((row) => [...row])
  let sign = 1
  let previous = Polynomial.one
  for (let k = 0; k < n - 1; k++) {
    // A nonzero pivot, swapping rows if need be; none means the determinant is 0
    if (m[k]![k]!.isZero()) {
      const swap = m.findIndex((row, i) => i > k && !row[k]!.isZero())
      if (swap < 0) return Polynomial.zero
      ;[m[k], m[swap]] = [m[swap]!, m[k]!]
      sign = -sign
    }
    const pivot = m[k]![k]!
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) {
        const numerator = m[i]![j]!.mul(pivot).sub(m[i]![k]!.mul(m[k]![j]!))
        const quotient = numerator.divideExact(previous)
        if (quotient === undefined) throw new Error('Bareiss elimination: inexact division')
        m[i]![j] = quotient
      }
      m[i]![k] = Polynomial.zero
    }
    previous = pivot
  }
  const det = m[n - 1]![n - 1]!
  return sign < 0 ? det.neg() : det
}

/** The matrix without row `row` and column `column`. */
export function minorMatrix(matrix: PolynomialMatrix, row: number, column: number): Polynomial[][] {
  return matrix.filter((_, i) => i !== row).map((r) => r.filter((_, j) => j !== column))
}

/** The (row, column) cofactor: the signed determinant of the minor. */
export function cofactor(matrix: PolynomialMatrix, row: number, column: number): Polynomial {
  const det = determinant(minorMatrix(matrix, row, column))
  return (row + column) % 2 === 0 ? det : det.neg()
}
