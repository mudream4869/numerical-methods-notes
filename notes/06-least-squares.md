# 6. Continuous Least Squares

## 6.1 Polynomial

Given $f \in C[a, b]$, find a polynomial $P(x) = \sum_{i=0}^{n} a_i x^i$ minimizing

$$
E(a_0, \dots, a_n) = \int_a^b \big( f(x) - P(x) \big)^2 \, dx.
$$

Setting the partial derivatives to zero:

$$
\frac{\partial E}{\partial a_i} =
-2 \int_a^b x^i f(x) \, dx + 2 \sum_{j=0}^{n} a_j \int_a^b x^{i+j} \, dx = 0,
\qquad i = 0, \dots, n.
$$

This is a system of $n + 1$ linear equations (the normal equations).

## 6.2 Linear Independence and Orthogonality

Use a linearly independent set of functions $\lbrace \phi_i(x) \rbrace$ with
$\int_a^b \phi_i(x) \phi_j(x) \, dx = 0$ for $i \neq j$, and let $P(x) = \sum_i a_i \phi_i(x)$.
Then $\frac{\partial E}{\partial a_i} = 0$ decouples into

$$
a_i \int_a^b \phi_i(x)^2 \, dx = \int_a^b f(x) \phi_i(x) \, dx.
$$

## 6.3 Trigonometric Polynomial Approximation

Use $\lbrace \phi_i(x) \rbrace_{i = 0, \dots, 2n-1}$:

$$
\begin{aligned}
\phi_0(x)     &= \tfrac{1}{2} \\
\phi_k(x)     &= \cos kx, \quad k = 1, \dots, n \\
\phi_{n+k}(x) &= \sin kx, \quad k = 1, \dots, n - 1
\end{aligned}
$$

They are linearly independent and orthogonal on $[-\pi, \pi]$,
so we can approximate $f \in C[-\pi, \pi]$ by

$$
S_n(x) = \frac{a_0}{2} + a_n \cos nx + \sum_{k=1}^{n-1} \big( a_k \cos kx + b_k \sin kx \big).
$$

## 6.4 Discrete Analog of Fourier Series

Integrals can be hard to compute, so use a discrete set of points instead:
$\lbrace (x_j, y_j) \rbrace_{j = 0, \dots, 2m-1}$ with $x_j = -\pi + \frac{j}{m} \pi$.

Minimize

$$
E = \sum_{j=0}^{2m-1} \big( y_j - S_n(x_j) \big)^2,
$$

where (for $n < m$) $S_n$ takes the slightly different form

$$
S_n(x) = \frac{a_0 + a_n \cos nx}{2} + \sum_{k=1}^{n-1} \big( a_k \cos kx + b_k \sin kx \big).
$$

The solution is

$$
\begin{aligned}
a_k &= \frac{1}{m} \sum_{j=0}^{2m-1} y_j \cos k x_j, \quad k = 0, \dots, n \\
b_k &= \frac{1}{m} \sum_{j=0}^{2m-1} y_j \sin k x_j, \quad k = 1, \dots, n - 1
\end{aligned}
$$
