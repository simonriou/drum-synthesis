# Problem Formulation

Let the domain of interest be the cube with dimensions $L_x \times L_y \times L_z = 1 \times 1 \times 1$. Let $h$ and $k$ be the spatial and temporal resolutions, such that $1 = N_x h = N_y h = N_z h$. We assume the membrane is rigidly clamped at its boundaries, so Neumann boundary conditions are used. The drum cavity is also assumed rigid.

# Finite Difference Methods

## Membrane Scheme

Let $w$ be the membrane displacement, $w = w(x, y, t) = w_{l, m}^n$ where $x = hl$, $y = hm$, and $t = kn$. The scheme is:
$$
\delta_{tt} w = \gamma^2 \delta_{\Delta_{2}} w - 2 \sigma_{0} \delta_{t} w
$$
where:
- $\gamma = c / L$ ($c$ is wave speed, $L = 1$ here)
- $\sigma_{0} = 6 \log 10 / T_{60}$ ($T_{60}$ is the $60$ dB decay time in seconds)

Expanding, the recursion is:
$$
w_{l, m}^{n+1} = \frac{2}{2k\sigma_{0} + 1} w_{l, m}^n + \frac{k^2\gamma^2}{2k\sigma_{0} + 1} \underline{\mathcal{D}}_{2} w_{l, m}^n + \frac{2\sigma_{0}k - 1}{2\sigma_{0}k + 1} w_{l, m}^{n-1}
$$

For stability, require $\lambda \leq 1/\sqrt{2}$ with $\lambda = \gamma k / h$ (Courant number).

$\underline{\mathcal{D}}_{2}$ is the 2D Laplacian matrix under Neumann conditions, block-structured as follows (block size $N_y+1 \times N_y+1$): diagonal blocks (edge or interior) are:
$$
\begin{pmatrix}
    -2 & 1 \\
    1 & -3 & \ddots \\
    & \ddots & \ddots & \ddots\\
    & & \ddots & \ddots & 1\\
    & & & 1 & -3 & 1 \\
    & & & & 1 & -2
    \end{pmatrix}
    \quad \text{or} \quad
    \begin{pmatrix}
    -3 & 1 \\
    1 & -4 & \ddots \\
    & \ddots & \ddots & \ddots\\
    & & \ddots & \ddots & 1\\
    & & & 1 & -4 & 1 \\
    & & & & 1 & -3
\end{pmatrix}
$$

## Cavity Scheme

Let $\Psi$ be the acoustic field variation inside the drum cavity:
$$
\Psi = \Psi(x, y, z, t) = \Psi_{l, m, p}^n
$$
with $x = lh$, $y = mh$, $z = ph$, $t = kn$.

Assume sound propagation follows the 3D wave equation, approximated as:
$$
\delta_{tt} \Psi = \gamma^2 \delta_{\Delta_{3}}\Psi
$$

### Stability Condition

Using the ansatz $\Psi_{l, m, p}^{n} = z^{n}\exp \left[ j h ( l\beta_{x} + m \beta_{y} + p\beta_{z} ) \right]$, the stability condition is derived as:
$$
\lambda = \gamma k / h \leq \frac{1}{\sqrt{3}}
$$

### Finite Difference Scheme

Expanding gives the recursion:
$$
\Psi_{l, m, p}^{n+1} = 2 \Psi_{l, m, p}^n + \lambda^2 \underline{\mathcal{D}}_{3} \Psi_{l, m, p}^{n} - \Psi_{l, m, p}^{n-1}
$$
$\underline{\mathcal{D}}_{3}$ is the 3D Laplacian matrix under Neumann conditions, structured using $\underline{\mathcal{D}}_{2}$ blocks. Super/sub-diagonals are identities; diagonal blocks are $\underline{\mathcal{D}}_{2} - I$ (interior) or $\underline{\mathcal{D}}_{2} - 2I$ (first/last).

# Modal Analysis

We compute PDE solutions using modal analysis.

## Membrane

Assume a solution to the lossless 2D wave equation:
$$
u(x, y, t) = e^{ j\omega t } U(x, y)
$$
Plugging into the PDE:
$$
-\omega^2 U(x, y) = \gamma^2 \nabla^2 U(x, y)
$$
Assume $U(x, y) = X(x) Y(y)$, with Neumann boundary conditions on $[0,1]^2$:
$$
\frac{X''(x)}{X(x)} + \frac{Y''(y)}{Y(y)} = -k^2 \implies \lambda^2 + \mu^2 = k^2
$$
Solve:
$$
\begin{cases}
X''(x) + \lambda^2 X(x) = 0, \quad X'(0) = X'(1) = 0 \\
Y''(y) + \mu^2 Y(y) = 0, \quad Y'(0) = Y'(1) = 0
\end{cases}
$$
General solution:
$$
X(x) = \cos(n\pi x), \quad Y(y) = \cos(m\pi y)
$$
for $n, m \in \mathbb{N}$.

Thus,
$$
k_{n, m}^2 = \pi^2(n^2 + m^2), \quad \omega_{n, m} = \gamma \pi \sqrt{ n^2 + m^2 }
$$
The general solution:
$$
U(x, y) = \sum_{n=0}^{\infty} \sum_{m=0}^{\infty} A_{n, m} \cos(n\pi x)\cos(m\pi y)
$$
Coefficients $A_{n, m}$ (for $U(x, y, 0) = f(x, y)$):
$$
A_{n, m} = 4 \int_{0}^{1} \int_{0}^{1} f(x, y) \cos(n\pi x) \cos(m\pi y) \, dx \, dy
$$

Now, with damping:
$$
\frac{ \partial^2 u }{ \partial t^2 } = \gamma^2 \nabla^2 u - 2 \sigma_{0} \frac{ \partial u }{ \partial t }
$$
Assume
$$
u(x, y, t) = \sum_{n=0}^{\infty} \sum_{m=0}^{\infty} A_{n, m}(t) \Phi_{n, m}(x, y)
$$
with $\Phi_{n, m}(x, y) = \cos(n\pi x)\cos(m\pi y)$.

For each mode $(n, m)$:
$$
A_{n, m}''(t) + 2\sigma_{0}A_{n, m}'(t) + \gamma^2 k_{n, m}^2 A_{n, m}(t) = 0
$$
This is a damped harmonic oscillator:
$$
s^2 + 2\sigma_{0}s + \gamma^2 k_{n, m}^2 = 0 \implies s = -\sigma_{0} \pm \sqrt{ \sigma_{0}^2 - \gamma^2k_{n, m}^2 }
$$
- Under-damped: $\gamma^2k_{n, m}^2 > \sigma_{0}^2$ (oscillatory decay)
- Critically damped: $\gamma^2k_{n, m}^2 = \sigma_{0}^2$
- Over-damped: $\gamma^2k_{n, m}^2 < \sigma_{0}^2$

Assuming under-damped:
$$
A_{n, m}(t) = A_{n, m} e^{ -\sigma_{0} t } \cos(\omega_{n, m} t), \quad \omega_{n, m} = \sqrt{ \gamma^2k_{n, m}^2 - \sigma_{0}^2 }
$$
with $A_{n, m}$ as above.

Final solution:
$$
\boxed{
u(x, y, t) = \sum_{n=0}^{\infty} \sum_{m=0}^{\infty} A_{n, m} e^{ -\sigma_{0} t } \cos(\omega_{n, m} t) \cos(n\pi x) \cos(m\pi y)
}
$$
