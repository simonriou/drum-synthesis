# Problem formulation
Let the domain of interest be the cube of dimensions $Lx \times Ly \times Lz = 1 \times 1 \times 1$. Let $h$ and $k$ be the spatial and temporal resolutions, such that $1 = Nx \times h = Ny \times h = Nz \times h$. We assume that the membrane is rigidly clamped at its boundaries, meaning Neumann boundary conditions will be used. Similarly, the drum cavity is assumed to be rigid too. 
# Finite difference methods
## Membrane scheme
Let $w$ be the displacement on the membrane, $w = w(x, y, t) = w_{l, m}^n$ where $x = hl$, $y = hm$ and ${} t = kn {}$. The scheme we are using is
$$\begin{align*}
\delta_{tt} w = \gamma^2 \delta_{\Delta_{2}} w - 2 \sigma_{0} \delta_{t.} w
\end{align*}$$
where
- $\gamma = c / L$ with $c$ the celerity of the wave and $L$ a characteristic length. In our case, $L = 1$.
- $\sigma_{0} = 6 \log 10 / T_{60}$ with $T_{60}$ the $60$ dB decay time (seconds)
After expanding the scheme, we arrive at the recursion
$$\begin{align*}
w_{l, m}^{n+1} = \frac{2}{2k\sigma_{0} + 1} u_{l, m}^n + \frac{k^2\gamma^2}{2k\sigma_{0} + 1} \underline{\mathcal{D}}_{2} w_{l, m}^n + \frac{2\sigma_{0}k - 1}{2\sigma_{0}k + 1} w_{l, m}^{n-1}
\end{align*}$$
It is important to note that, in order for the stability of this scheme to be valid, one must have $\lambda \leq 1 / \sqrt{ 2 }$ with $\lambda = \gamma k / h$ the Courant number.
Here, $\underline{\mathcal{D}}_{2}$ refers to the matrix operator of the 2D Laplacian approximation under Neumann conditions. It can be decomposed as blocks, and is of shape $Nx + 1 \times Nx + 1$ blocks, each of size $Ny + 1 \times Ny + 1$. Still thinking of it as blocks, the super and sub diagonals are identity blocks. The diagonal is composed of blocks of the following shape:
$$
\begin{pmatrix}
-2 & 1 \\
1 & -3 & \ddots \\
& \ddots & \ddots & \ddots\\
& & \ddots & \ddots & 1\\
& & & 1 & -3 & 1 \\
& & & & 1 & -2
\end{pmatrix} \ \text{ or } \begin{pmatrix}
-3 & 1 \\
1 & -4 & \ddots \\
& \ddots & \ddots & \ddots\\
& & \ddots & \ddots & 1\\
& & & 1 & -4 & 1 \\
& & & & 1 & -3
\end{pmatrix}
$$
depending if they correspond to edge or interior points.
## Cavity scheme
In the context of the 3D study, we will denote $\Psi$ the variation in the acoustic field inside of the drum cavity. We will also write
$$\begin{align*}
\Psi = \Psi(x, y, z, t) = \Psi_{l, m, p}^n
\end{align*}$$
where $x = lh, y = mh, z = ph$ and $t = kn$.
The drum cavity being a relatively small volume, we assume that the propagation of sound in it follows the 3D wave equation, which we approximate using the following scheme:
$$\begin{align*}
\delta_{tt} \Psi = \gamma^2 \delta_{\Delta_{3}}\Psi
\end{align*}$$
### Stability condition
To determine the 3D stability condition, we consider the ansatz $\Psi_{l, m, p}^{n} = z^{n}\exp \left[ j h \left( l\beta_{x} + m \beta_{y} + p\beta_{z} \right) \right]$. The left-hand term of the equation becomes
$$\begin{align*}
\delta_{tt} \Psi &= \frac{1}{k^2} \left[ z^{n+1} - 2z^{n} + z^{n-1} \right] \exp \left[ j h \left( l\beta_{x} + m \beta_{y} + p\beta_{z} \right) \right]\\
&= \frac{z^{n}}{k^2}\left[ z - 2 + z^{-1} \right] \exp \left( j \Phi \right) 
\end{align*}$$
The right-hand term is
$$\begin{align*}
\delta_{\Delta_{3}} \Psi &= \left( \delta_{xx} + \delta_{yy} + \delta_{zz} \right)\Psi \\
&= z^n \left( A_{x} + A_{y} + A_{z} \right)
\end{align*}$$
Developing the $x$-wise term gives
$$\begin{align*}
A_{x} &= \frac{1}{h^2} \left[ e^{ jh\left( l+1 \right) \beta_{x} } - 2 e^{ jhl\beta_{x} } + e^{ jh\left( l-1 \right) \beta_{x} } \right] \\
&= \frac{e^{ jhl\beta_{x} }}{h^2}\left( 2 \cos(h \beta_{x}) - 2 \right)
\end{align*}$$
From which we can deduce that
$$\begin{align*}
&A_{y} =\frac{e^{ jhm\beta_{y} }}{h^2}\left( 2 \cos(h \beta_{y}) - 2 \right) & A_{z} = \frac{e^{ jhp\beta_{z} }}{h^2}\left( 2 \cos(h \beta_{z}) - 2 \right)
\end{align*}$$
Meaning that
$$\begin{align*}
\delta_{\Delta_{3}} \Psi = \frac{2}{h^2}\left[ \cos(h\beta_{x}) + \cos(h\beta_{y}) + \cos(h\beta_{z}) -3 \right] 
\end{align*}$$
And so the PDE becomes
$$\begin{align*}
\frac{z^{n}}{k^2}\left[ z - 2 + z^{-1} \right]  = \frac{2\gamma^2}{h^2}\left[ \cos(h\beta_{x}) + \cos(h\beta_{y}) + \cos(h\beta_{z}) -3 \right] z^{n}
\end{align*}$$
Now denoting $\lambda = \gamma k / h$, we can write this equation as
$$\begin{align*}
&z + z^{-1} - 2 = 2 \lambda^2 \left[ \underbrace{\cos(h\beta_{x}) + \cos(h\beta_{y}) + \cos(h\beta_{z})}_{\zeta} - 3 \right] \\
\implies & z + z^{-1} = 2 + 2\lambda^2\left( \zeta - 3 \right)
\end{align*}$$
To ensure the stability of the scheme, we need $|z| < 1$, i.e. $z + z^{-1} \in [-2, 2]$, i.e. $2 + 2\lambda^2\left( \zeta - 3 \right) \in [-2, 2]$. The worst-case scenario happens for $\zeta = -3$ and we then have
$$\begin{align*}
&2 -12 \lambda^2 \geq -2\\
\implies &\lambda^2 \leq \frac{1}{3}\\
\implies &\boxed{\lambda \leq \frac{1}{\sqrt{ 3 }}}
\end{align*}$$
### Finite difference scheme
Just as in the 2D case, expanding the scheme leads to the following recursion
$$\begin{align*}
\Psi_{l, m, p}^{n+1} = 2 \Psi_{l, m, p}^n + \lambda^2 \underline{\mathcal{D}}_{3} \Psi_{l, m, p}^{n} - \Psi_{l, m, p}^{n-1}
\end{align*}$$
where, similarly as before, $\underline{\mathcal{D}}_{3}$ refers to the matrix operator of the 3D Laplacian approximation under Neumann conditions. This matrix can be expressed using $\underline{\mathcal{D}}_{2}$, as it is composed of $Nz + 1 \times Nz + 1$ blocks, each of size $(Ny + 1)(Nx + 1) \times (Ny + 1)(Nx + 1)$. The super and sub-diagonals of $\underline{\mathcal{D}}_{3}$ are identities, and the diagonal blocks are variations of $\underline{\mathcal{D}}_{2}$: interior blocks are 2D Laplacian matrix from which we subtracted the identity, while the first and last blocks are $\underline{\mathcal{D}}_{2} - 2I$.
This is however only valid if we consider the cavity alone. As soon as we introduce coupling between the membrane and the cavity, the matrix operator of the 3D Laplacian must be updated to fit the new boundary conditions.
## Coupling
The air pressure just beneath the membrane exerts a force on it. Explicitly, we have that
$$\begin{align*}
\delta_{t.} w_{l, m}^{n} = -\lim_{ p \to N_{z} } \delta_{z.} \Psi_{l, m, p}^{n}
\end{align*}$$
In other words, the velocity of the membrane is equal to the opposite of the vertical velocity of the acoustic field evaluated right below the membrane, i.e. at the topmost layer of the 3D cavity. Differentiating this relation leads to
$$\begin{align*}
\frac{1}{2k}\left( w_{l, m}^{n+1} - w_{l, m}^{n-1} \right) = -\frac{1}{h}\left( \Psi_{l, m, N_{z}}^n - \Psi_{l, m, N_{z} - 1}^{n}\right) 
\end{align*}$$
This relation will be used to determine the effect of the vibration of the membrane on the acoustic field inside the cavity. 
We now need a way to determine the inverse process, i.e. what is the influence of the acoustic field inside the drum on the displacement of the membrane. To do so, we use the following coupling term
$$\begin{align*}
f_{l, m}^{-, n} = \rho \lim_{ p \to N_{z} } \delta_{t.} \Psi_{l, m, p}^n
\end{align*}$$
which can be written as
$$\begin{align*}
f_{l, m}^{-, n} = \frac{1}{2k}\left( \Psi_{., ., N_{z}}^{n+1} - \Psi_{., ., N_{z}}^{n-1} \right) 
\end{align*}$$

If we had to detail the process of computing $w^2$ from $w^{1}$, it would be
1. We know $w_{l, m}^1$ as it is the excitation (typically a raised cosine distribution). We consider that $w_{l, m}^0 = 0$. We consider that the acoustic field inside the cavity $\Psi$ is zero everywhere at $n = 0$.
2. The displacement of the membrane at time steps $1$ and $0$ allows us to determine the value of $\Psi_{., ., N_{z}}^1$ or $\Psi_{., ., N_{z}}^0$, I'm not sure (coupling conditions).
3. The 3D wave equation allows us to then determine $\Psi_{l, m, p}^2$ or $\Psi_{l, m, p}^1$, I don't know, for all $p$. The influence of the membrane is contained in the previous iteration of the acoustic field (Cf. step 2).
4. $\Psi_{l, m, p}^2$ and $\Psi_{l, m, p}^1$ allows us to determine $f^-(n = 1)$ the coupling term.
5. $f^-(n=1), w_{l, m}^0$ and $w_{l, m}^1$ are now used to determine $w_{l, m}^2$.
6. We repeat.
# Modal analysis
In this section, we are computing solutions to the PDEs using modal analysis.
## Membrane
Assuming a solution to the lossless 2D wave equation to be of the form
$$\begin{align*}
u(x, y, t) = e^{ j\omega t } U(x, y)
\end{align*}$$
One may write, plugging this solution back into the PDE, that
$$\begin{align*}
-\omega^2 U(x, y) = \gamma^2 \ \nabla^2 U(x, y)
\end{align*}$$
where ${} \nabla^2 {}$ refers to the 2D Laplacian operator. We will assume that $U$ is separable, meaning we can write it as $U(x, y) = X(x) \times Y(y)$, as well as the fact that $U(x, y)$ is defined over a rectangular domain with Neumann boundary conditions. We can re-write the previous equation as
$$\begin{align*}
&\nabla^2 U + k^2 U = 0 & k^2 = \frac{\gamma^2}{\omega^2}
\end{align*}$$
or, replacing $U$ by $X \times Y$,
$$\begin{align*}
\frac{X''(x)}{X(x)} + \frac{Y''(y)}{Y(y)} = -k^2 \implies \lambda^2 + \mu^2 = k^2
\end{align*}$$
This leads to a system of two equations to solve:
$$\begin{align*}
\begin{cases}
X''(x) + \lambda^2 X(x) = 0 & X'(0) = X'(1) = 0 & (1) \\
Y''(y) + \mu^2 Y(y) = 0 & Y'(0) = Y'(1) = 0 & (2)\\
\end{cases}
\end{align*}$$
We will detail the computations for $(1)$. The general form of the solution is
$$\begin{align*}
X(x) = A\cos(\lambda x) + B \sin(\lambda x)
\end{align*}$$
The boundary conditions give us that $X'(0) = 0 \implies B = 0$ and so that $X(x) = A \cos(\lambda x)$, and that $X'(1) = 0 \implies \lambda_{n} = n\pi$. The choice of $A$ is arbitrary, thus we will consider $A = 1$. The solutions are then
$$\begin{align*}
&X(x) = \cos(n\pi x) & Y(y) = \cos(m\pi y)
\end{align*}$$
This means that for every pair of integers ${} (n, m) \in \mathbb{N}^{*2} {}$, we have
$$\begin{align*}
&\lambda_{n} = n\pi & \mu_{m} = m\pi && k_{n, m}^2 = \lambda_{n}^2 + \mu_{m}^2 = \pi^2(n^2 + m^2)
\end{align*}$$
The general solution is thus
$$\begin{align*}
U(x, y) &= X(x) Y(y)\\
&= \sum_{n = 0}^{\infty} \ \sum_{m = 0}^{\infty} \ A_{n, m} \cos(n\pi x)\cos(m\pi y)
\end{align*}$$
And this solution satisfies the 2D wave equation is each of the modes $(n, m)$ satisfies
$$\begin{align*}
&\pi^2(n^2 + m^2) = \frac{\omega^2}{\gamma^2} &\text{or} &&k_{n, m}^2 = k^2
\end{align*}$$
This means that the pulsations of the system are
$$\begin{align*}
\omega_{n, m} = \gamma \pi \sqrt{ n^2 + m^2 }
\end{align*}$$
The coefficients $A_{n, m}$ depend on the initial state of the system. Assuming the general case where $U(x, y, 0) = f(x, y)$, then these coefficients are given by Fourier analysis:
$$\begin{align*}
A_{n, m} = 4 \int_{0}^{1} \int_{0}^{1} f(x, y) \cos(n\pi x) \cos(m\pi y) \ dx \ dy
\end{align*}$$
And the final solution is
$$\begin{align*}
|u(x, y, t)| = U(x, y)
\end{align*}$$
Now, we are interested in the case where there is a damping term to the equation:
$$\begin{align*}
\frac{ \partial^2 u }{ \partial t^2 } = \gamma^2 \nabla^2 u - 2 \sigma_{0} \frac{ \partial u }{ \partial t }
\end{align*}$$
From what we have done above, we can assume that
$$\begin{align*}
u(x, y, t) = \sum_{n = 0}^{\infty} \ \sum_{m = 0}^{\infty } \ A_{n, m}(t) \Phi_{n, m}(x, y)
\end{align*}$$
Where $\Phi_{n, m}(x, y) = \cos(n\pi x)\cos(m\pi y)$ which satisfies the Neumann boundary conditions over $[0, 1] \times [0, 1]$. Noticing that
$$\begin{align*}
&\nabla^2 \Phi_{n, m} = -k_{n, m}^2 \Phi_{n, m} & k_{n, m}^2 = \pi^2(n^2 + m^2)
\end{align*}$$
we can plug back the solution into the PDE, which means that for each mode $(n, m)$,
$$\begin{align*}
&A_{n, m}''(t) \Phi_{n, m} + 2\sigma_{0}A_{n, m}'(t) \Phi_{n, m} = -\gamma^2 k_{n, m}^2 \Phi_{n, m}A_{n, m}(t)\\
\implies &A_{n, m}''(t) + 2\sigma_{0}A_{n, m}'(t) + \gamma^2 k_{n, m}^2 A_{n, m}(t) = 0
\end{align*}$$
Which is the form of a classic damped harmonic oscillator. Its characteristic equation is
$$\begin{align*}
s^2 + 2\sigma_{0}s + \gamma^2 k_{n, m}^2 = 0 \implies s = -\sigma_{0} \pm \sqrt{ \sigma_{0}^2 - \gamma^2k_{n, m}^2 }
\end{align*}$$
And from there, three cases arise:
- Under-damped: $\gamma^2k_{n, m}^2 > \sigma_{0}^2$ and the behaviour is oscillatory with decay
- Critically damped: $\gamma^2k_{n, m}^2 = \sigma_{0}^2$
- Over-damped: $\gamma^2k_{n, m}^2 < \sigma_{0}^2$

We will assume the first case and thus can write
$$\begin{align*}
&A_{n, m} = C_{n, m} e^{ -\sigma_{0} t } \cos(\omega_{n, m} t) & \omega_{n, m} = \sqrt{ \gamma^2k_{n, m}^2 - \sigma_{0}^2 }
\end{align*}$$
Assuming zero initial velocity (as the initial condition will most likely be a raised cosine distribution), we have that $A_{n, m}'(0) = 0$ and thus
$$\begin{align*}
u(x, y, 0) = f(x, y) = \sum_{n = 0}^{\infty} \ \sum_{m = 0}^{\infty} \ A_{n, m}(0) \Phi_{n, m}
\end{align*}$$
From which we can deduce that $C_{n, m} = A_{n, m}(0)$, and thus that
$$\begin{align*}
A_{n, m}(t) &= e^{ -\sigma_{0}t } \cos(\omega_{n, m}t) \times 4 \int_{0}^{1} \int_{0}^{1} f(x, y) \cos(n\pi x) \cos(m\pi y) \ dx \ dy \\
&= e^{ -\sigma_{0} t} \cos(\omega _{n, m} t) A_{n, m}
\end{align*}$$
Which gives us the final solution
$$\begin{align*}
\boxed{u(x, y, t) = \sum_{n = 0}^{\infty} \ \sum_{m = 0}^{\infty } \ A_{n, m} e^{ -\sigma_{0} t} \cos(\omega_{n, m}t) \cos(n\pi x) \cos(m\pi y)}
\end{align*}$$
