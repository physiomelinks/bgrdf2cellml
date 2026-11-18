/**
 * Exact rational numbers on BigInt, always in lowest terms with a positive
 * denominator, so equal values are equal objects field by field.
 */

export interface Rational {
  readonly n: bigint
  readonly d: bigint
}

const gcd = (a: bigint, b: bigint): bigint => {
  a = a < 0n ? -a : a
  b = b < 0n ? -b : b
  while (b !== 0n) [a, b] = [b, a % b]
  return a
}

/** `n/d` in lowest terms. */
export function rational(n: bigint | number, d: bigint | number = 1n): Rational {
  let num = BigInt(n)
  let den = BigInt(d)
  if (den === 0n) throw new RangeError('Rational with a zero denominator')
  if (den < 0n) {
    num = -num
    den = -den
  }
  const g = gcd(num, den)
  return g > 1n ? { n: num / g, d: den / g } : { n: num, d: den }
}

export const ZERO: Rational = { n: 0n, d: 1n }
export const ONE: Rational = { n: 1n, d: 1n }

/**
 * A decimal as written (`"2.5e-3"`, `"0.726"`) exactly, not by way of a
 * binary float: 0.1 is 1/10.
 */
export function fromDecimal(text: string): Rational {
  const match = /^([+-]?)(\d*)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$/.exec(text.trim())
  if (match === null || (match[2] === '' && (match[3] ?? '') === '')) throw new SyntaxError(`Not a decimal number: ${text}`)
  const [, sign, whole, fraction = '', exponent = '0'] = match
  let n = BigInt(`${whole}${fraction}` || '0')
  let d = 10n ** BigInt(fraction.length)
  const e = Number(exponent)
  if (e > 0) n *= 10n ** BigInt(e)
  else if (e < 0) d *= 10n ** BigInt(-e)
  return rational(sign === '-' ? -n : n, d)
}

/** A JavaScript number exactly, through its shortest decimal form. */
export function fromNumber(value: number): Rational {
  if (!Number.isFinite(value)) throw new RangeError(`Not a finite number: ${value}`)
  if (Number.isInteger(value) && Math.abs(value) <= Number.MAX_SAFE_INTEGER) return rational(value)
  return fromDecimal(String(value))
}

export const isZero = (a: Rational): boolean => a.n === 0n
export const isOne = (a: Rational): boolean => a.n === 1n && a.d === 1n
export const equal = (a: Rational, b: Rational): boolean => a.n === b.n && a.d === b.d
export const sign = (a: Rational): -1 | 0 | 1 => (a.n < 0n ? -1 : a.n > 0n ? 1 : 0)

export const add = (a: Rational, b: Rational): Rational => (a.d === b.d ? rational(a.n + b.n, a.d) : rational(a.n * b.d + b.n * a.d, a.d * b.d))
export const neg = (a: Rational): Rational => ({ n: -a.n, d: a.d })
export const sub = (a: Rational, b: Rational): Rational => add(a, neg(b))
export const mul = (a: Rational, b: Rational): Rational => rational(a.n * b.n, a.d * b.d)

export function div(a: Rational, b: Rational): Rational {
  if (b.n === 0n) throw new RangeError('Division by zero')
  return rational(a.n * b.d, a.d * b.n)
}

/** `a` to an integer power. */
export function pow(a: Rational, k: number): Rational {
  if (!Number.isInteger(k)) throw new RangeError(`Not an integer power: ${k}`)
  if (k < 0) return pow(div(ONE, a), -k)
  const e = BigInt(k)
  return { n: a.n ** e, d: a.d ** e }
}

export function compare(a: Rational, b: Rational): number {
  const x = a.n * b.d
  const y = b.n * a.d
  return x < y ? -1 : x > y ? 1 : 0
}

export const toNumber = (a: Rational): number => (a.d === 1n ? Number(a.n) : Number(a.n) / Number(a.d))

export const toString = (a: Rational): string => (a.d === 1n ? `${a.n}` : `${a.n}/${a.d}`)

/** The integer if `a` is one, else `undefined`. */
export const asInteger = (a: Rational): number | undefined => (a.d === 1n ? Number(a.n) : undefined)

/** The greatest common divisor of two rationals: gcd of numerators over lcm of denominators, positive. */
export function gcdRational(a: Rational, b: Rational): Rational {
  if (a.n === 0n) return { n: b.n < 0n ? -b.n : b.n, d: b.d }
  if (b.n === 0n) return { n: a.n < 0n ? -a.n : a.n, d: a.d }
  const dg = gcd(a.d, b.d)
  return rational(gcd(a.n, b.n), (a.d / dg) * b.d)
}
