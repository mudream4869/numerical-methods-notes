# 4. Nonlinear Equations

## 4.1 Nonlinear Equations of One Variable

Given a continuous function of one variable, find a root of it.

### Bisection Method

If $f(a) f(b) < 0$, there is a root in $(a, b)$.
Let $m = \frac{a + b}{2}$; check the sign of $f(a) f(m)$ to decide whether the root lies in
$(a, m)$ or $(m, b)$. Repeat to shrink the interval.

### Newton's Method

Approximate $f$ by its tangent line at $x_i$ and solve for its zero:

$$
\begin{aligned}
&f(x_i) + f'(x_i)(x - x_i) = 0 \\
\Rightarrow\ &x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}
\end{aligned}
$$

Note that it may not converge.

## 4.2 Convergence Rate

We show that, near a root, Newton's method converges faster than bisection.
Let $r$ be a root of $f$, i.e. $f(r) = 0$.

> **Definition 4.1.** If $\lim_{k \to \infty} \frac{|x_{k+1} - r|}{|x_k - r|} = \delta < 1$,
> the sequence has **linear convergence**.

Bisection converges linearly (the error bound halves each step).

> **Definition 4.2.** If $\lim_{k \to \infty} \frac{|x_{k+1} - r|}{|x_k - r|^2} = \delta < \infty$,
> the sequence has **quadratic convergence**.

### Newton's Method: Quadratic Convergence

Assume $f$ is continuously differentiable and $f'$ is Lipschitz continuous with constant $\alpha$:
$|f'(x) - f'(y)| \le \alpha |x - y|$.

> **Theorem 4.1.** If the Newton iterates $\lbrace x_k \rbrace \to x^\ast$ and $f'(x^\ast) \neq 0$, then
>
> 1. $f(x^\ast) = 0$;
> 2. with $\beta = |f'(x^\ast)|^{-1}$, $\exists L$ s.t. $\forall k \ge L$,
>    $|x_{k+1} - x^\ast| \le \alpha \beta \, |x_k - x^\ast|^2$.

*Proof.* Define $e(y, x) = f(y) - f(x) - f'(x)(y - x)$. Then

$$
\begin{aligned}
|e(y, x)| &= \left| \int_x^y \big( f'(t) - f'(x) \big) \, dt \right| \\
&\le \left| \int_x^y \alpha |t - x| \, dt \right| \\
&= \frac{\alpha}{2} |y - x|^2.
\end{aligned}
$$

(1) By the iteration, $f(x_k) = f'(x_k)(x_k - x_{k+1})$.
Since $f'$ is continuous and $x_k - x_{k+1} \to 0$, $f(x^\ast) = \lim_{k \to \infty} f(x_k) = 0$.

(2) Since $f'$ is continuous and $f'(x^\ast) \neq 0$, $\exists L$ s.t. $\forall k \ge L$,
$|f'(x_k)| \ge \frac{1}{2} |f'(x^\ast)|$, i.e. $|f'(x_k)^{-1}| \le 2\beta$.

From $0 = f(x^\ast) = f(x_k) + f'(x_k)(x^\ast - x_k) + e(x^\ast, x_k)$:

$$
\begin{aligned}
|x_{k+1} - x^\ast|
&= \left| x_k - x^\ast - \frac{f(x_k)}{f'(x_k)} \right| \\
&= \left| x_k - x^\ast - \frac{-f'(x_k)(x^\ast - x_k) - e(x^\ast, x_k)}{f'(x_k)} \right| \\
&= |f'(x_k)^{-1}| \, |e(x^\ast, x_k)| \\
&\le 2\beta \cdot \frac{\alpha}{2} |x_k - x^\ast|^2 = \alpha \beta \, |x_k - x^\ast|^2. \qquad \blacksquare
\end{aligned}
$$

## 4.3 Newton's Method for Multiple Variables and Equations

Consider $m$ equations $f: \mathbb{R}^n \to \mathbb{R}^m$ in $n$ variables.

When $m = n$, Newton's method generalizes directly with the Jacobian $J$:

$$
x_{k+1} = x_k - J(x_k)^{-1} f(x_k)
$$

When $m \neq n$, $J$ is $m \times n$, not square, hence not invertible.

Instead minimize $g(x) = \Vert f(x) \Vert^2$, i.e. solve $\nabla g(x) = 2 J(x)^T f(x) = 0$.
This is $n$ equations in $n$ variables, and its Jacobian $\nabla^2 g$ is $n \times n$, so Newton's method applies.
(Approximating $\nabla^2 g \approx 2 J^T J$ gives the Gauss-Newton method.)

## 4.4 Fixed-Point Iteration

For example, the system

$$
\begin{aligned}
f_1(x_1, x_2) &= x_1^3 + 10 x_1 - x_2 - 5 = 0 \\
f_2(x_1, x_2) &= x_2^3 - 10 x_2 + x_1 + 1 = 0
\end{aligned}
$$

can be rewritten as $x = g(x)$:

$$
\begin{aligned}
x_1 &= g_1(x_1, x_2) = -0.1 x_1^3 + 0.1 x_2 + 0.5 \\
x_2 &= g_2(x_1, x_2) = 0.1 x_2^3 + 0.1 x_1 + 0.1
\end{aligned}
$$

> **Theorem 4.2.** If $p = g(p)$, $g$ is continuously differentiable,
> $\Vert J_g(p) \Vert < 1$ for some induced norm, and $x_0$ is sufficiently close to $p$,
> then $x_{k+1} = g(x_k)$ converges to $p$.
