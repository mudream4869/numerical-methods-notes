# 3. Norms and Iterative Methods

## 3.1 Norms of Vectors and Matrices

A **norm** satisfies:

$$
\begin{aligned}
\Vert A \Vert &\ge 0, \text{ and } \Vert A \Vert = 0 \iff A = 0 \\
\Vert A + B \Vert &\le \Vert A \Vert + \Vert B \Vert \\
\Vert \alpha A \Vert &= |\alpha| \, \Vert A \Vert
\end{aligned}
$$

Some vector norms:

1. 2-norm on $\mathbb{R}^2$: $\Vert (x, y) \Vert = \sqrt{x^2 + y^2}$
2. $p$-norm on $\mathbb{R}^n$: $\Vert x \Vert_p = \left( \sum_{i=1}^{n} |x_i|^p \right)^{1/p}$
3. max-norm on $\mathbb{R}^n$: $\Vert x \Vert_\infty = \max_i |x_i|$

Note that when $p < 1$, the $p$-norm is not a norm (triangle inequality fails).

We could treat a matrix as a "long vector", but that loses matrix-specific properties.
A matrix norm is usually the **induced norm**

$$
\Vert A \Vert = \max_{\Vert v \Vert = 1} \Vert A v \Vert,
$$

which satisfies $\Vert Av \Vert \le \Vert A \Vert \Vert v \Vert$ and is useful for linear systems.

### Perturbing the RHS

Consider $Ax = b$ and perturb the RHS: $A(x + \delta x) = b + \delta b$, so $A \, \delta x = \delta b$.

$$
\begin{aligned}
\Vert \delta x \Vert = \Vert A^{-1} \delta b \Vert &\le \Vert A^{-1} \Vert \Vert \delta b \Vert \\
\Vert b \Vert = \Vert A x \Vert &\le \Vert A \Vert \Vert x \Vert
\end{aligned}
$$

Combining:

$$
\frac{\Vert \delta x \Vert}{\Vert x \Vert} \le \Vert A \Vert \Vert A^{-1} \Vert \frac{\Vert \delta b \Vert}{\Vert b \Vert}
$$

### Perturbing the Matrix

Perturb $A$: $(A + \delta A)(x + \delta x) = b$.

$$
\begin{aligned}
&\delta x = -A^{-1} \, \delta A \, (x + \delta x) \\
\Rightarrow\ &\Vert \delta x \Vert \le \Vert A^{-1} \Vert \Vert \delta A \Vert \Vert x + \delta x \Vert \\
\Rightarrow\ &\frac{\Vert \delta x \Vert}{\Vert x + \delta x \Vert} \le \Vert A \Vert \Vert A^{-1} \Vert \frac{\Vert \delta A \Vert}{\Vert A \Vert}
\end{aligned}
$$

$\kappa(A) = \Vert A \Vert \Vert A^{-1} \Vert$ is called the **condition number** of $A$.

## 3.2 Iterative Methods for Linear Systems

For $Ax = b$ in $\mathbb{R}^n$, Gaussian elimination takes $O(n^3)$ time.
An iterative method improves the current approximate solution step by step.
We want (number of iterations) $\times$ (cost per iteration) $\ll O(n^3)$.

### Jacobi and Gauss-Seidel Methods

Let $D = \operatorname{diag}(A)$. The **Jacobi method** iterates

$$
D \, x^{(k+1)} = (D - A) \, x^{(k)} + b,
$$

i.e. $x_i^{(k+1)} = \frac{1}{A_{ii}} \left( b_i - \sum_{j \neq i} A_{ij} x_j^{(k)} \right)$.
Each iteration costs $O(n^2)$ ($O(n)$ per component). For $n = 3$:

$$
\begin{aligned}
x_1^{(k+1)} &= \frac{b_1 - A_{12} x_2^{(k)} - A_{13} x_3^{(k)}}{A_{11}} \\
x_2^{(k+1)} &= \frac{b_2 - A_{21} x_1^{(k)} - A_{23} x_3^{(k)}}{A_{22}} \\
x_3^{(k+1)} &= \frac{b_3 - A_{31} x_1^{(k)} - A_{32} x_2^{(k)}}{A_{33}}
\end{aligned}
$$

**Gauss-Seidel** uses the newly updated components as soon as they are available:

$$
\begin{aligned}
x_1^{(k+1)} &= \frac{b_1 - A_{12} x_2^{(k)} - A_{13} x_3^{(k)}}{A_{11}} \\
x_2^{(k+1)} &= \frac{b_2 - A_{21} x_1^{(k+1)} - A_{23} x_3^{(k)}}{A_{22}} \\
x_3^{(k+1)} &= \frac{b_3 - A_{31} x_1^{(k+1)} - A_{32} x_2^{(k+1)}}{A_{33}}
\end{aligned}
$$

Note that the iteration may diverge (both converge if $A$ is strictly diagonally dominant).

### Conjugate Gradient Method

Assume $A$ is symmetric positive definite. Solving $Ax = b$ is equivalent to

$$
\min_x f(x) = \frac{1}{2} x^T A x - b^T x,
$$

since $\nabla f(x) = Ax - b = -r$ (the residual).

**Steepest descent** moves along $-\nabla f(x_k) = r_k$:
$x_{k+1} = x_k + \alpha r_k$. Minimizing over the step size $\alpha$ is a one-variable problem:

$$
\alpha = \frac{r_k^T r_k}{r_k^T A r_k}.
$$

The denominator is positive since $A$ is positive definite ($x^T A x > 0$ for $x \neq 0$).
This algorithm can be **very slow**.

Instead, use search directions $p_1, \dots, p_n$ with $x_n \in x_0 + \operatorname{span}\lbrace p_1, \dots, p_n \rbrace$.
Writing $\sum_{i < n} y_i p_i = P_{n-1} y$, the problem becomes

$$
\min_{y, \alpha} f(x_0 + P_{n-1} y + \alpha p_n).
$$

If $p_n$ is **$A$-conjugate** to the previous directions ($p_n^T A p_i = 0$ for $i < n$,
i.e. $p_n \perp \operatorname{span}\lbrace A p_1, \dots, A p_{n-1} \rbrace$), the cross term vanishes and the problem splits:

$$
f(x_0 + P_{n-1} y + \alpha p_n) = f(x_0 + P_{n-1} y) - \alpha \, p_n^T r_0 + \frac{\alpha^2}{2} p_n^T A p_n.
$$

The $y$ part is the previous step's problem (solved by $x_{n-1}$), so

$$
x_n = x_{n-1} + \alpha_n p_n, \qquad \alpha_n = \frac{p_n^T r_0}{p_n^T A p_n}.
$$

The conjugacy constraint still leaves freedom, so choose $p_n$ closest to the residual $r_{n-1} = b - A x_{n-1}$:

$$
p_n = \arg\min_{p \,\perp\, \operatorname{span}\lbrace A p_1, \dots, A p_{n-1} \rbrace} \Vert p - r_{n-1} \Vert
$$

TODO: show that this gives

$$
p_n \in \operatorname{span}\lbrace r_{n-1}, p_{n-1} \rbrace.
$$
