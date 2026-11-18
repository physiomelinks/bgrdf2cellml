# @physiomelinks/steady-state

Derives steady-state flux expressions for enzyme and transporter cycles, symbolically, from their mass-action or bond-graph equations: the reductions of Hunter, Ai & Nickerson (2025) for GLUT2 and SGLT1, done automatically.

```js
import { deriveSteadyState, verify } from '@physiomelinks/steady-state'

const spec = { label: 'GLUT2', cycleStates: ['q_1', 'q_2', 'q_3', 'q_4'], flux: 'v_R2', fast: ['kappa_R3', 'kappa_R4'] }
const result = deriveSteadyState(equations, spec) // equations: math-ast Equation[]
result.common.flux        // v = κ/Φ · {q̄_Ao − q̄_Ai}, in named lumped parameters
result.common.definitions // kappa_m_GLUT2, k_Ao_GLUT2, …, Phi_GLUT2 in the model's parameters
verify(equations, spec, result).maxRelativeError // a numeric check against the full cycle
```

It is built from these parts:

- `Q`: exact rationals on BigInt, reading decimals such as `0.726` exactly;
- `Polynomial`: sparse multivariate polynomials with rational coefficients (add, multiply, substitute, evaluate, content, exact division);
- `determinant` and `cofactor`: Bareiss fraction-free elimination over polynomials, which gives a cycle's King–Altman state weights as cofactors of its rate matrix.
- `normalise`: turns a model's expressions into exact quotients of polynomials. It inlines definitions, cancels `RT·ln(K·q)/RT` so that `exp(u_f/RT)` becomes the mass-action product (`K₁q₁·(K₂q₂)²`), and writes every other exponent term, such as `z·F·u_m/RT`, as an *exp atom*. Atoms whose exponents are multiples of each other share one atom, so `exp(F·u/(2RT))² / exp(F·u/RT)` is 1. Constructs it can't represent (piecewise, other functions, a logarithm outside an exponent) are reported as a `SteadyStateError` with an issue code.
- `fractionToAst` and `polynomialToAst` write the results back as `@physiomelinks/math-ast` trees, and `evaluateAst` evaluates a tree numerically, for checking.
- `solveCycle`: the steady-state flux of a linear cycle. Given the model's equations, the cycle's states and the flux to report, it:
  - builds the rate matrix;
  - checks that the states are conserved and enter linearly;
  - takes each state's King–Altman weight as a cofactor of the matrix;
  - returns the flux as an exact quotient of polynomials, in the parameters, the held species, the total `q_tot` and the exp atoms.

  When it can't reduce a model it returns issues (`no-conservation`, `nonlinear-in-states`, `disconnected-cycle`, `unknown-variable`) instead.
- **Rapid equilibrium:** give `solveCycle` the rate constants of fast reactions (`fast: ['kappa_R3', 'kappa_R4']`) and it takes their limit by Cha's method.
  - States joined by fast reactions are pooled in their equilibrium proportions, and the pooled cycle is solved.
  - A flux through a fast reaction is rewritten in terms of slow ones using the steady-state balances.
  - This reproduces GLUT2's eqs 19–20 (Hunter et al. 2025) and Kir's eq 19 (sodium homeostasis paper).
- `deriveSteadyState` writes the result in the **common transporter form** of the sodium homeostasis paper (eq 24):

  `v = κ/Φ · {R · Q_drive/Q_drive^eq − Q_eq · P · exp(zFu/RT)}`

  - The bracket is a *difference* (`q̄_i − q̄_o`, with `Φ = 1 + Σ m/k_m`) or a *ratio* (`(q̄_Na,o/q̄_Na,i)²·q̄_Glc,o/q̄_Glc,i − e^{2Fu/RT}`).
  - It can be written in concentrations (`concentrations: { q_Ao: 'V_o' }`).
  - An optional driver reaction (ATP hydrolysis) gets its own ratio.
  - It reports the net charge `z` and the reversal potential.
  - `relations` adds constraints such as Kir's charge balance before solving.
  - A flux with more than one cycle (slippage) is still solved exactly, but is reported as not having this form.
- `verify` checks a result numerically, independently of the algebra: at random values it solves the model's own cycle equations with floating-point linear algebra.
