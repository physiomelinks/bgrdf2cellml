# @physiomelinks/steady-state

Derives steady-state flux expressions for enzyme and transporter cycles, symbolically, from their mass-action or bond-graph equations: the reductions of Hunter, Ai & Nickerson (2025) for GLUT2 and SGLT1, done automatically.

Work in progress. This first part is the exact algebra it is built on:

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
