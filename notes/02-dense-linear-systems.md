# 2. Solving Dense Linear Systems

## 2.1 LU Factorization

Goal: solve the linear system $Ax = b$.

### Triangular Systems

A triangular matrix is a square matrix whose entries above (or below) the diagonal are zero.
For example,

$$
\begin{bmatrix}
1 & 2 & 3 \\
0 & 4 & 5 \\
0 & 0 & 6
\end{bmatrix}
$$

is an upper triangular matrix.

When $A$ is triangular, $Ax = b$ is solved by forward/back substitution in
$O(n^2)$ time, where $A$ is $n \times n$.

Observe that

$$
\begin{bmatrix} 3 & 5 \\ 6 & 7 \end{bmatrix}
=
\begin{bmatrix} 1 & 0 \\ 2 & 1 \end{bmatrix}
\begin{bmatrix} 3 & 5 \\ 0 & -3 \end{bmatrix},
$$

so we can solve

$$
\begin{bmatrix} 3 & 5 \\ 6 & 7 \end{bmatrix} x = \begin{bmatrix} 9 \\ 4 \end{bmatrix}
$$

by splitting $A$ into two triangular matrices and solving two triangular systems:

$$
\begin{bmatrix} 1 & 0 \\ 2 & 1 \end{bmatrix}
\begin{bmatrix} 3 & 5 \\ 0 & -3 \end{bmatrix} x = \begin{bmatrix} 9 \\ 4 \end{bmatrix}.
$$

### Gaussian Elimination

Gaussian elimination can be written as $U = M_{n-1} M_{n-2} \cdots M_1 A$,
where $U$ is upper triangular and each $M_i$ is a (unit lower triangular) row-operation matrix.
Hence $L = M_1^{-1} \cdots M_{n-1}^{-1}$ is also lower triangular, and $A = LU$.

### Row Interchange (Pivoting)

Not every nonsingular $A$ can be factorized by plain Gaussian elimination. For instance,

$$
\begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}
$$

has a zero pivot. Interchanging rows fixes this: $PA = LU$ for some permutation matrix $P$.

TODO: partial pivoting details.

## 2.2 Cholesky Factorization

> **Theorem 2.1.** A symmetric matrix $A$ is positive definite $\Leftrightarrow$
> $A = LL^T$ for some lower triangular $L$ with positive diagonal entries.

*Proof.*

($\Rightarrow$) Induction on $n$. For $n = 1$, $A = [\alpha]$ with $\alpha > 0$, so $L = [\sqrt{\alpha}]$.

For $n > 1$, write

$$
A = \begin{bmatrix} \alpha & v^T \\ v & B \end{bmatrix},
$$

where $\alpha = e_1^T A e_1 > 0$. Let $S = B - \frac{v v^T}{\alpha}$ (Schur complement).
$S$ is positive definite: for $y \neq 0$, take $x = \begin{bmatrix} -v^T y / \alpha \\ y \end{bmatrix}$, then $x^T A x = y^T S y > 0$.
By induction $S = L_S L_S^T$. Then

$$
A =
\begin{bmatrix} \sqrt{\alpha} & 0 \\ \frac{v}{\sqrt{\alpha}} & L_S \end{bmatrix}
\begin{bmatrix} \sqrt{\alpha} & \frac{v^T}{\sqrt{\alpha}} \\ 0 & L_S^T \end{bmatrix},
$$

since the lower-right block of the product is $\frac{v v^T}{\alpha} + L_S L_S^T = B$.

($\Leftarrow$) $L$ is nonsingular, so for $x \neq 0$, $x^T A x = \Vert L^T x \Vert^2 > 0$. $\blacksquare$
