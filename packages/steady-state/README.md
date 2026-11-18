# @physiomelinks/steady-state

Derives steady-state flux expressions for enzyme and transporter cycles, symbolically, from their mass-action or bond-graph equations: the reductions of Hunter, Ai & Nickerson (2025) for GLUT2 and SGLT1, done automatically.

Work in progress. This first part is the exact algebra it is built on:

- `Q`: exact rationals on BigInt, reading decimals such as `0.726` exactly;
- `Polynomial`: sparse multivariate polynomials with rational coefficients (add, multiply, substitute, evaluate, content, exact division);
- `determinant` and `cofactor`: Bareiss fraction-free elimination over polynomials, which gives a cycle's King–Altman state weights as cofactors of its rate matrix.
