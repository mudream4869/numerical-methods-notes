# 5. Interpolation and Approximation

## 5.1 Fitting a Function

Given $(x_1, f(x_1)), \dots, (x_n, f(x_n))$ with distinct $x_i$, find a function passing through all of them.

### Lagrange Interpolation

$$
P(x) = \sum_{i=1}^{n} f(x_i) \, L_i(x), \qquad
L_i(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}
$$

Since $L_i(x_j) = 1$ if $i = j$ and $0$ otherwise, $P(x_i) = f(x_i)$. $P$ has degree at most $n - 1$.

### Spline

The Lagrange polynomial has too high a degree for many points (it oscillates).
Instead, use a different cubic polynomial on each interval $[x_i, x_{i+1}]$,
requiring the pieces and their first and second derivatives to be continuous at the interior nodes.

## 5.2 Function Approximation

Given data $(x_i, y_i)$ that need not be fitted exactly, find $P$ minimizing

$$
\sum_i \big( y_i - P(x_i) \big)^2.
$$

TODO
