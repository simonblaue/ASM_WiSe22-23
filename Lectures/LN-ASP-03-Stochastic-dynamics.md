---
title: "Lecture Notes on Advanced Statistical Physics 03 "
header-includes:
- \usepackage{braket}

output:
  pdf_document: 
       latex_engine: xelatex

---

[TOC]


## 3. Stochastic dynamics and temporal coarse-graining

Motivation:
: How do we include interaction with heat bath in dynamics? 

![Colloid ~1µm in solvent]()

### 3.1 Caldeira-Leggett model, Langevin equation

Guess of Hamiltonian:
$$
\begin{align*}
 H = \frac{p^2}{2m} + V(x) + \sum_{i} \left(\frac{p_i^2}{2m} + \frac{m_i\omega_i^2}{2}(x_i-x)^2\right)
\end{align*}
$$

Can rewrite this as

$$
\begin{align*}
    H = \frac{p^2}{2m} + V_1(x) + \sum_{i} \left(\frac{p_i^2}{2m} + \frac{m_i\omega_i^2}{2} x_i^2\right) +\underbrace{\sum_i\left( \frac{m_i\omega_i^2}{2} x^2  \right)}_{V_2(x)} - \underbrace{\sum_i\left( {m_i\omega_i^2} x_ix  \right)}_{\text{bilinear coupling}},
\end{align*}$$

with ${V'}(x) = V_1 + V_2$

Aim: Eliminate ('coarse-grain) the $x_i, p_i$ by solving their dynamics explicitly

$$
	\dot p = =  -\frac{\partial H}{\partial x} = -V'(x)+ \sum_i m_i\omega_i^2(x_i-x) = m \ddot x 
$$

$$
\dot p_i =  -\frac{\partial H}{\partial x_i} = - m_i\omega_i^2(x_i-x) = m_i \ddot x_i
$$

Solution of equations for $x_i(t)$:

$$
\begin{align*}
    x_i(t) = x_i(0) \cos(\omega_i t)+ \sin(\omega_i t) \frac{p_i(0)}{m_i\omega_i} + \int_0^t dt' \omega_i\sin(\omega_i(t-t'))x(t')
\end{align*}$$

Check:

$
    \dot x_i(t) = - \omega_i x_i(0) \sin(\omega_i t)+ \cos(\omega_i t) \frac{p_i(0)}{m_i} + \int_0^t dt' \omega_i^2\cos(\omega_i(t-t'))x(t')
$

$
    \ddot x_i(t) = - \omega_i^2x_i(0) \cos(\omega_i t)- \sin(\omega_i t) \frac{p_i(0) \omega_i}{m_i} - \int_0^t dt' \omega_i^2\sin(\omega_i(t-t'))x(t') + \omega_i^2 x(t) = -\omega_i^2 x_i(t)+ \omega_i^2 x(t)
$

✓

Initial conditions:

$x_i(0) = x_i(0)\cdot 1 + 0 + 0$ ✓
$\dot x_i(0) = 0 + \frac{p_i(0)}{m_i} + 0$ ✓

Integral can be integrated by parts:

$\int_0^t {\sin(\omega_i(t-t'))}x(t') = \underbrace{\frac{x(t')}{\omega_i}\cos(\omega_i(t-t'))|_0^t}_{=\frac{1}{\omega_i}(x(t)-x(0)\cos(\omega_i(t-t')))} - \int_0^tdt'\frac{1}{\omega_i}\cos(\omega_i(t-t'))\dot x(t')$

Insert $x_i(t)$ int equation for $x(t)$:

$$
\begin{align*}
    m\ddot x = -V'(x) + \sum_im_i\omega_i^2(x_i-x) \\
    = -V'(x) + \sum_i m_i \omega_i^2 \underbrace{[(x_i(0)-x(0))\cos(\omega_i t) + \frac{p_i(0)}{m_i\omega_i}\sin(\omega_i t) + x(t)-x(t)]}_{\xi(t)} \\
    - \int_0^t dt' \underbrace{\sum_i m_i\omega_i^2 \cos(\omega_i(t-t'))}_{\gamma(t-t')}\dot x(t')
\end{align*}$$

Result ==Generalized Langevin equation==:

$$
\begin{align*}
    m\ddot x = -V'(x) - \int_0^t dt' \underbrace{\gamma(t-t')}_{\text{friction kernal}} \dot x(t') + \underbrace{\xi(t)}_{\text{effective noise/random force}}
\end{align*}$$

To rewrite the friction kernel, introduce the "spectral density" of oszillators

$ J(\omega) = \sum_i \frac{\pi}{2} m_i \omega_i^3 \delta(\omega-\omega_i) $:

$$\begin{align*}
\gamma(t) = \sum_i m_i\omega_i^2\cos(\omega_i t) = \frac{2}{\pi} \sum_i \frac{\pi}{2}m_i\omega_i^3 \frac{1}{\omega_i}\cos(\omega_i t) \\
= \frac{2}{\pi} \sum_i \frac{\pi}{2} m_i \omega_i^3 \int_0^{\infin} \frac{d\omega}{\omega} \cos(\omega t) \delta(\omega-\omega_i) = \frac{2}{\pi} \int_0^{\infin}d\omega \frac{J(\omega)}{\omega}\cos(\omega t) \\
\frac{2}{\pi}\frac{1}{2}\int_{-\infin}^{\infin} \frac{d\omega}{|\omega|} \frac{J(|\omega|)}{|\omega|}\cos(\omega t) = \frac{2}{\pi}\frac{1}{4}\int_{-\infin}^{\infin} \frac{d\omega}{|\omega|} J(|\omega|)(e^{i\omega t}+ e^{-i\omega t}) \\
\int_{-\infin}^{\infin} \frac{d\omega}{2\pi} \frac{J(|\omega|)}{|\omega|}(e^{i\omega t}+ e^{-i\omega t})
\end{align*}$$

Consider $J(\omega)$ now as a smooth function (many oscillators) "Ohmic bath".

$J(\omega)=\gamma_0 \omega$ up to some $\Omega$ or with soft cutoff:

$J(\omega)=\frac{\gamma_0\omega}{1+(\omega/\Omega)^2}$ gives,

$\gamma(t) = \gamma_0 \Omega e^{-\Omega|t|}$ as frictional kernel

---

__TODO VL. 14.11__

---

### 3.2 Fokker-Planck equation

- So far: heat bath generates friction + noise
- Generelized Langevin equation
  - $m\ddot x + \gamma \dot x = -V'(x) + \xi$
  - overdamped Langevin: $\gamma \dot x = -V'(X)+ \xi$
- Aim: Translate into equation for $P(x,t), A(x,t)$

> How does $x$ change in some $\Delta t$?  Ex. $V=\frac{k}{2} x^2$

$\gamma \dot x = -kx + \xi \qquad \qquad$ , $\braket{\xi(t)\xi(t')}=2T\gamma \delta(t-t')$

$$\begin{align*}
\Delta x = \frac{-k}{\gamma} \int_t^{t+\Delta t} x(t') dt' + \frac{1}{\gamma} \xi_{\text{tot}}\\
\simeq \frac{-k}{\gamma}x(t)\Delta t +\frac{1}{\gamma}\xi_{\text{tot}}
\end{align*}$$

Reminder: $\braket{\xi^2_{\text{tot}}} = 2T\gamma \Delta t$
So $x'=x(t) +\Delta x$ has a Gaussian distribution:
$$\begin{align*}
    P(x'|x,\Delta t) = \mathscr{N}(x'|\underbrace{x-\frac{k}{\gamma}x\Delta t}_{\text{mean}}, \underbrace{\frac{2T}{\gamma}\Delta t}_{\text{variance}})
\end{align*}$$

Start now from generic $P(x'|x,\Delta t)$ Then

Forward Kolmogorov equation:
:   $$\begin{align*}
    P(x, t+\Delta t) = \int dx' P(x|x',\Delta t)P(x',t)
    \end{align*}$$

key input here is lack of memory/ __Markovianity__

> How do observables evolove, say $A(x,t)$, $x$ initial condition ?

Need to define $A$ now as ==average== value of observable given initial condition $x$ (and time $t$)

Evoulution of $A(x,t)$:

$$\begin{align*}
    A(x,t+\Delta t) = \int dx' A(x',t) P(x'|x,\Delta t)
\end{align*}$$

Convert into differentiable equation using "Kramers-Moyal" expansion: By Taylor,
$$\begin{align*}
    A(x,t+\Delta t) = \int dx' [A(x,t)+\sum_{n\geq1}\frac{(x'-x)^n}{n!}\frac{\partial^nA(x,t)}{\partial x^n}]P(x'|x,\Delta t) \\
    = A(x,t) + \sum_{n\geq 1} \frac{\partial^n A(x,t)}{n! \partial x^n} \int dx' (x'-x)^n P(x'|x,\Delta t)
\end{align*}$$
So
$$
    \frac{A(x,t+\Delta t)- A(x,t)}{\Delta t} = \sum_{n\geq 1}\frac{1}{n!}\frac{\partial^n A(x,t)}{\partial x^n} \underbrace{\int dx' \frac{(x'-x)^n}{\Delta t} P(x'|x,\Delta t)}_{\text{jump moment} M^{(n)}(x)}
\end{align*}$$

Take $\Delta t\rightarrow 0$, define
$$\begin{align*}
D^{(n)}(x)= \lim_{\Delta t\rightarrow 0}\frac{1}{n!}M^{(n)}(x)
\end{align*}$$

$$\begin{align*}
    \frac{\partial}{\partial t}A(x,t) = \sum_{n\geq 1}D^{(n)}(x)\frac{\partial^n}{\partial x^n}A(x,t)
\end{align*}$$

Can now convert to equation for $P(x,t)$ Use
$$\begin{align*}
\braket{A} = \int dx A(x,t)P(x,0) = \int dx A(x,0)P(x,t)
\end{align*}$$

Take $\frac{\partial}{\partial t}$ at $t=0$

$$\begin{align*}
    \int dx \frac{\partial A}{\partial t}(x,t)P(x,0) = \int A(x,0) \frac{\partial}{\partial t}P(x,t) \\
    \Leftrightarrow \int dx \sum_{n\geq 1} D^{(n)}(x)\left(\frac{\partial^nA}{\partial x^n}(x,t)\right)P(x,0) = \int dx A(x,0) \frac{\partial}{\partial t}P(x,t) \\
    \Leftrightarrow \int dx \sum_{n\geq 1} A(x,0) \left((-1)^n\frac{\partial^n}{\partial x^n}\right)(D^{(n)}(x)P(x,0)) = \int dx A(x,0) \frac{\partial}{\partial t}P(x,0)
\end{align*}$$
Therby (works for arbitrary $t$),
$$\begin{align*}
    \frac{\partial P}{\partial t}(x,t) = \sum_{n\geq 1} (-1)^n \frac{\partial^n}{\partial x^n}(D^{(n)}(x)P(x,t))
\end{align*}$$

#### 3.2.1 Harmonic oszillator example

$$\begin{align*}
    M^{(1)}(x) = \frac{\braket{x'-x}}{\Delta t}= -\frac{k}{\gamma} \frac{x \Delta t}{\Delta t} = -\frac{kx}{\gamma}, D^{(1)}(x) = M^{(1)}(x) = -\frac{kx}{\gamma}
\end{align*}$$

$$\begin{align*}
    M^{(2)}(x) = \frac{\braket{(x'-x)^2}}{\Delta t} = \frac{1}{\Delta t}\left (\underbrace{\left(-\frac{k}{\gamma}x\Delta t\right)^2}_{\text{mean}^2} + \underbrace{\frac{2T}{\gamma} \Delta T}_{\text{variance}} \right) = \frac{2T}{\gamma} + \mathscr{O}(\Delta t)
\end{align*}$$
So
$$\begin{align*}
    D^{(2)}(x) = \lim_{\Delta t\rightarrow 0} \frac{1}{2!}  M^{(2)}(x) = \frac{T}{\gamma}
\end{align*}$$
All higher $D^{(n)}(x)$ vanish. Intuitive argument:

$$\begin{align*}
    x'-x = \underbrace{\mathscr{O}(\Delta t)}_{mean} +   \underbrace{\mathscr{O}((\Delta t)^{1/2})}_{std. dev}
\end{align*}$$
So $(x'-x)^n\sim\mathscr{O}((\Delta t)^{n/2} )$, so $M^{(n)}(x) = \mathscr{O}((\Delta t)^{\frac{n}{2}-1})$ and $M^{(n)}(x)\rightarrow 0$ for $\Delta t \rightarrow 0$ if $n\geq 3$. Hence:

$$\begin{align*}
    \frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}\left(-\frac{kx}{\gamma}
    P\right)+\frac{\partial^2}{\partial x^2}\left(\frac{T}{\gamma}P\right)
\end{align*}$$

#### 3.2.2 Derivation (Kramers-Moyal) extends to higher dimensions

$$\begin{align*}
    \frac{\partial}{\partial t}P(\vec{x},t) = \sum_{n\geq 1}\sum_{i_1...i_n}(-\frac{\partial}{\partial x_{i_1}})(-\frac{\partial}{\partial x_{i_2}})...(-\frac{\partial}{\partial x_{i_n}})(D^{(n)}_{i_1...i_n}(\vec{x})P(\vec{x},t))
\end{align*}$$

For Langevin dynamics: $\dot{\vec x}= \vec v(\vec x)+ \vec \xi$ need only $n=1,2$

$$\begin{align*}
    \frac{\partial}{\partial t} P(\vec x,t) = -\sum_{i}\frac{\partial}{\partial x_i}(v_i(\vec x)P(\vec x,t)) + \sum_{i,j}\frac{\partial^2}{\partial x_i\partial x_j}(D_{ij}(\vec x)P(\vec x,t))
\end{align*}$$

Called the ==Fokker-Planck== equation, reduces to Liouville for $D_{ij}=0$

If $D_{i,j}(\vec x)$ is constant, noise is __additive__, otherwise "multiplicative" noise
Noise then has to be interpreted carefully, derivation so far uses "Itô" interpretation. (noise depends on $\vec x(t)$, not $\vec x(t+\Delta t)$)

__Interprete special cases of FP:__

Overdamped dynamics of particles $\vec v_i$, potential $V(\vec r)$, FP becomes ==Smoluchowski equation==.

$$\begin{align*}
    \gamma \dot{\vec r_i} = -\vec \nabla_i V+ \xi_i \quad \text{or} \quad \gamma \dot{\vec r_{i\alpha}} = -\frac{\partial V}{\partial r_{i\alpha}} + \xi_{i\alpha}
\end{align*}$$

with $\braket{\xi_{i\alpha}(t)\xi_{i\beta}(t')} = 2T\gamma \delta_{ij}\delta_{\alpha\beta}\delta(t-t')$.

$$\begin{align*}
    \frac{\partial P}{\partial t} = -\sum_{i\alpha}\frac{\partial}{\partial r_{i\alpha}}(-\frac{1}{\gamma}\frac{\partial V}{\partial r_{i\alpha}}P) + \sum_{i\alpha, j\beta} \frac{\partial}{\partial r_{i\alpha}} \frac{\partial}{\partial r_{j\beta}}(\frac{T}{\gamma}\delta_{ij}\delta{\alpha\beta}P)
\end{align*}$$

$$\begin{align*}
    = +\sum_i \vec \nabla_i \cdot (\frac{1}{\gamma}(\vec\nabla_i V)P) + \frac{T}{\gamma}\sum_i(\vec\nabla_i \cdot \vec\nabla_i)P \\
    = \frac{1}{\gamma} \vec \nabla (\vec \nabla V) P + \frac{T}{\gamma} \nabla^2 P, \vec \nabla = (\vec \nabla_1, ..., \vec \nabla_N)
\end{align*}$$

---

#### Dynamics with inertial in 1D

$$\begin{align*}
    m\ddot x = -V'(x)-\gamma\dot x+\xi, \qquad \braket{\xi(t)\xi(t')}=2T\gamma\delta(t-t')
\end{align*}$$
or

$$\begin{align*}
    m\dot x=p \\
    \dot p = -V'(x)-\frac{\gamma}{m}p+\xi
\end{align*}$$

So setting $\vec x =(x,p)$: $v_1 = \frac{p}{m}$, $v_2=-V'(x)-\frac{\gamma}{m}p$, $\xi_1=0$, $\xi_2=\xi$ so $D_{11}=D_{12}=0$, $D_{22}=\frac{1}{2!}2T\gamma = T\gamma$ gives Klein-Kramer equation:

$$\begin{align*}
    \partial_t P(x,p,t)= -\partial_x(\frac{p}{m}P)-\partial_p((-V'(x)-\frac{\gamma}{m}p)P)+\partial_p^2(T\gamma P)
\end{align*}$$

---

### 3.3 Probability currents, stationary states, detailed balanced

> Fokker-Planck can be wirtten as continuity equation !

$$\begin{align*}
    \partial_t P = -\sum_i \partial_{x_i} j_i, \quad j_i = v_iP-\sum_j\partial_{x_j}(D_{ij}P)
\end{align*}$$

Can interpret $j_i/P$ as velocity

$$\begin{align*}
    \frac{j_i(\vec x,t)}{P}=\underbrace{v_i(\vec{x})}_{\text{deterministic drift}} - \underbrace{\sum_j \partial_{x_j}D_{ij}}_{\text{spurious drift}} - \underbrace{\sum_j D_{ij}\underbrace{\frac{1}{P} \partial_{x_j}P}_{\partial_{x_j} \ln(P)}}_{\text{velocity of diffusive current}}
\end{align*}$$

#### 3.3.1 Stationary state of Fokker-Planck: $\partial_t P^{\text{ss}}=0$

$$\begin{align*}
    \sum_i \partial_{x_i}j_i^{\text{ss}} = 0,\quad \text{current is divergence free}
\end{align*}$$
Stationary states with $j_i^{\text{ss}}\equiv 0$ are called ==equilibrium==. A system with an equilibrium state is said to be ==microscopical reversible== or to obey ==detailed balance==

---

##### 3.3.1.0 Smoluchowski as example

$$\begin{align*}
 j = -\frac{1}{\gamma}\vec\nabla V-\frac{T}{\gamma}\vec\nabla P
\end{align*}$$

"Guess" Boltzman Distribution for equilibirum $P^{\text{ss}}(\vec{r})=\frac{1}{z}e^{-\beta V(\vec r)}$, gives

$$\begin{align*}
    j^{\text{ss}} = -\frac{1}{\gamma}(\vec \nabla V)\frac{1}{z}e^{-\beta V(\vec r)}-\frac{T}{\gamma} \underbrace{\vec\nabla \frac{1}{z}e^{-\beta V(\vec r)}}_{=\frac{1}{z}(-\beta \vec\nabla V)e^{-\beta V(\vec r)}} = 0
\end{align*}$$

> Note: Detailed Balanced (DB) only refers to stationary state
> System will "relax" into equilibrium state if started out of equilibrium.

E.g. for Smoluchowski, can show that free energy:

$$\begin{align*}
    F[P] = \int d\vec r P(\vec r) V(\vec r) + \underbrace{T \int d\vec r P(\vec r) \ln P(\vec r)}_{= -TS[P]}
\end{align*}$$

decreases monotonically in time. ($-d_t F \propto \int ... j^2$)

---

Stationary states with $j^{\text{ss}}\neq 0$ can ocure with e.g.

- external forcing
- explicit breaking of DB, __active matter__

Example with external forcing: Kramers problem
(escape from metastable potential well)

![[Pasted image 20221121132820.png]]
  

Average number of particles: $\bar N$
particle injection rate: $\Gamma = \bar N j$
$j =$ escape rate $= \frac{\Gamma}{\bar N}$

Alternatively, normalize to $\bar N= 1$, then $j=\Gamma $

Write FP in terms of $\ln P$:

$$\begin{align*}
	\gamma \partial_t\ln P = \frac{\gamma}{P}\partial_t P= \frac{1}{P} (\partial_x V'(x)e^{\ln P}) + T\frac{1}{P}(\partial_x^2 e^{\ln P}) \\
	= \frac{1}{P}(V''P+V'(\ln P)'P) + \frac{T}{P}(\underbrace{\partial_x(\ln P)' e^{\ln P}}_{[(\ln P)''+((\ln P)')^2]P}) \\ = V'' + V'(\ln P)' + T(\ln P)'' + T((\ln P)')^2
\end{align*}$$

Simplify for low $T$: Expect $P$ to scale as $e^{-U(x)/T}$.
So $\ln P = -U/T + \text{const.}$, leading terms $[0(\frac{1}{T})]$ in stationary state condition $0=V'' + ...$ are
$$\begin{align*}
    0 = \frac{1}{T}(V'(-U')+(-U')^2) = \frac{1}{T}U'(U'-V')
$$

So inside metastable potential well:

$$\begin{align*}
U' = V', U(x) = V(x)\text{const.} = V(x)-V(x_{\text{min}}) \quad (\text{for} P=e^{-U/T} )
\end{align*}$$

While outside $U'=0, U=\text{const.}$ so $j=\mathscr{O}(1) e^{-(V(x_\text{barrier})-V(x_\text{min}))/T}$, hence

$$\begin{align*}
    j\propto e^{-\Delta V/T}
\end{align*}$$

---

### 3.4 Correlation functions

So far: 

$$\begin{align*}
\partial_t P = -\sum_i \partial_{x_i}(v_i P) + \sum_{ij} \partial_{x_{ij}}^2(D_{ij}P)= L^\dagger P
\end{align*}$$

(dual):

$$\begin{align*}
\partial_t A = \sum_i v_i \partial_{x_i}(A) + \sum_{ij} D_{ij}\partial_{x_{ij}}^2(A)= L A
\end{align*}$$

> $L \neq -L^\dagger$

Correlation functions: 
 
$$\begin{align*}
    C_{AB}(t,t') = \int d\vec x A(\vec x,t)B(\vec x,t')P(\vec x, 0)
\end{align*}$$

but this misses correlations from stochastic dynamics. except of $t'=0$ (or $t=0$)
