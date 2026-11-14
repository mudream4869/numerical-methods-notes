# 1. Floating Point Operation

## 1.1 Normalized Representation

Consider a floating point system with base $\beta$, precision $p$ and
exponent range $[e_{\min}, e_{\max}] \cap \mathbb{Z}$.
A number in this system can be written in **normalized representation**

$$
\pm\, d_0.d_1 d_2 \dots d_{p-1} \times \beta^{e}, \qquad d_0 \neq 0,
$$

i.e. $a \times \beta^e$ with $1 \le |a| < \beta$. For example:

1. $0.1 = 1.00 \times 10^{-1}$ ($\beta = 10$, $p = 3$)
2. $0.1 \approx 1.1001 \times 2^{-4}$ ($\beta = 2$, $p = 5$)

This representation cannot represent $0$, since $|a| \ge 1$.
The natural representation for $0$ is $1.0 \times \beta^{e_{\min} - 1}$.

## 1.2 Relative Errors and Units in the Last Place (ulps)

The **ulp** (unit in the last place) is used as a unit of error.
For $d_0.d_1 \dots d_{p-1} \times \beta^e$, one ulp is $\beta^{-p+1} \cdot \beta^e$.

The largest error when representing a real number $r$:

$$
E(r) = \min |a \times \beta^e - r|
$$

Let $r = d_0.d_1 d_2 \dots d_{p-1} d_p \dots \times \beta^e$.
Choose $s = d_0.d_1 d_2 \dots d_{p-2} k \times \beta^e$, where

$$
k =
\begin{cases}
d_{p-1}     & \text{if } 0.d_p d_{p+1} \dots < \tfrac{1}{2} \\
d_{p-1} + 1 & \text{otherwise (with carry)}
\end{cases}
$$

Then $E(r) \le \frac{\beta}{2} \beta^{-p} \beta^e = \frac{1}{2}\,\text{ulp}$.

Since $|r| \ge \beta^e$, the relative error is bounded by

$$
\frac{E(r)}{|r|} \le \frac{E(r)}{\beta^e} \le \frac{\beta^{-p+1}}{2}.
$$

This bound is called the **machine epsilon** $\epsilon = \frac{\beta^{-p+1}}{2}$.

## 1.3 Guard Digits

To compute $x - y$, **compute then round** gives small error. Example ($p = 3$, $\beta = 10$):

$$
\begin{aligned}
x     &= 2.1500000000 \\
y     &= 0.0000000125 \\
x - y &= 2.1499999875
\end{aligned}
$$

which rounds to $2.15$.

**Round then compute** (shift $y$ to $p$ digits first) gets the same answer more cheaply:

$$
\begin{aligned}
x     &= 2.15 \\
y     &= 0.00 \\
x - y &= 2.15
\end{aligned}
$$

However, round then compute can cause a large relative error:

$$
\begin{aligned}
x      &= 1.00  \\
y      &= 0.999 \\
y'     &= 0.99  \quad (\text{truncated}) \\
x - y' &= 0.01  \\
x - y  &= 0.001
\end{aligned}
$$

Relative error $= \frac{|0.01 - 0.001|}{0.001} = 9$.

> **Theorem 1.1.** Using $p$ digits with base $\beta$ to compute $x - y$,
> the relative error can be as large as $\beta - 1$.

Large errors occur when $x$ and $y$ are close. A **single guard digit**
prevents this: *the adder/subtractor works with $p + 1$ digits.*

> **Theorem 1.2.** Using $p + 1$ digits (one guard digit) for $x - y$,
> the relative error is $< 2\epsilon$.

*Proof.* TODO. See Goldberg, *What Every Computer Scientist Should Know About
Floating-Point Arithmetic* (1991), Theorem 2.

## 1.4 Cancellation

Let $a = 1.22$, $b = 3.34$, $c = 2.28$ ($p = 3$). Exactly, $b^2 - 4ac = 0.0292$. But

$$
\begin{aligned}
b^2       &= 11.1556 \approx 11.2 \\
4ac       &= 11.1264 \approx 11.1 \\
b^2 - 4ac &\approx 0.1
\end{aligned}
$$

$b^2$ and $4ac$ already contain rounding errors, and subtracting two close
numbers exposes them. This is called **catastrophic cancellation**.

How to avoid it? Reformulate the expression.

### $\frac{-b + \sqrt{b^2 - 4ac}}{2a}$

When $b > 0$ and $b^2 \gg 4ac$, we have $\sqrt{b^2 - 4ac} \approx b$, so the numerator cancels.
Multiplying by the conjugate gives

$$
\frac{-b + \sqrt{b^2 - 4ac}}{2a} = \frac{2c}{-b - \sqrt{b^2 - 4ac}},
$$

which avoids this cancellation.
Note that not every cancellation can be removed (e.g. the one inside $b^2 - 4ac$).

### $x^2 - y^2$

When $x$ and $y$ are close, $(x - y)(x + y)$ is more accurate than $x^2 - y^2$.

## 1.5 Rounding

Round to nearest: e.g. $1.4 \to 1$, $1.6 \to 2$.
For ties, **round half to even**: $2.5 \to 2$, $3.5 \to 4$.

> **Theorem 1.3** (Reiser & Knuth). Let $x_0 = x$ and $x_n = (x_{n-1} - y) + y$.
> If $+$ and $-$ use round half to even, then either
> $x_n = x_0$ for all $n \ge 1$, or $x_n = x_1$ for all $n \ge 1$.

## 1.6 IEEE Standard

How to choose $\beta$ and $p$? Compare two systems using the same 4 bits:

1. $\beta = 16$, $p = 1$: $\epsilon = \frac{16}{2} \cdot 16^{-1} = \frac{1}{2}$
2. $\beta = 2$, $p = 4$: $\epsilon = \frac{2}{2} \cdot 2^{-4} = \frac{1}{16}$

A smaller base gives a smaller $\epsilon$ for the same storage; IEEE 754 uses $\beta = 2$.

## 1.7 Special Quantities

### Infinity

IEEE has $+0$, $-0$, $+\infty$, $-\infty$ to represent values too small or too
large for the format. Arithmetic follows limits, e.g.
$\frac{1}{+0} = +\infty$ since $\lim_{x \to 0^+} \frac{1}{x} = +\infty$.

### Not a Number

When an operation has no real result, the result is $\mathrm{NaN}$:

$$
\begin{aligned}
\sqrt{-4}       &= \mathrm{NaN} \\
\infty - \infty &= \mathrm{NaN}
\end{aligned}
$$

### $+0$ and $-0$

IEEE defines the comparison $+0 = -0$ to be true, which is convenient
(e.g. `if (x == 0)`). However the sign still matters:
$+0 = -0$ while $\frac{1}{+0} = +\infty \neq -\infty = \frac{1}{-0}$.

(In IEEE 754, $\log(\pm 0) = -\infty$ and $\log x = \mathrm{NaN}$ for $x < 0$.)

### Denormalized Numbers

To represent numbers smaller than $1.0 \times \beta^{e_{\min}}$:
$0.d_1 d_2 \dots d_{p-1} \times \beta^{e_{\min}}$ (gradual underflow).

## 1.8 Exceptions and Traps

Sometimes we want to handle special quantities ourselves instead of continuing.
The system should let the programmer choose: **use special values** or **trap**.

| Situation      | Result without trap         |
| -------------- | --------------------------- |
| overflow       | $\pm\infty$                 |
| underflow      | $0$ or a denormalized number |
| divide by zero | $\pm\infty$                 |
| invalid        | $\mathrm{NaN}$              |
| inexact        | $\mathrm{round}(x)$         |

### Example: Power Function

For $n < 0$, $x^n$ can be computed as `1 / pow(x, -n)` or `pow(1/x, -n)`.
Choose the one that does not trigger a spurious overflow or underflow.

### Example: arccos

TODO
