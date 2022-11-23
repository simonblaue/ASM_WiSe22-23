---
title: "Lecture Notes on Advanced Statistical Physics "
header-includes:
- \usepackage{braket}

output:
  pdf_document: 
       latex_engine: xelatex

---



## Expectations

- [ ] combinatroics
- [x] Foker-Plank eq.
- [x] stochastic processes, Langevin eq.
- [x] phase seperation/ transitions
- [x] non equilibrium physics
- [x] nonlinear problems
- [ ] reaction kinetics, network dynamics
- [x] information
- [x] quantum statistical physics, mostly classical
- [x] linear response theory, FDT
- [x] diagrammatic methods (Feynmann in statistical physics)
- [x] (field therory )

## Background

- thermodynamics eg. thermodynamic potentials, legendre transformations
- ensembles, ==Boltzman distribution==
- (detailed balance)
- density matrices
- (brownian motion)

---

# 1. Lecture Notes on Advanced Statistical Physics


- [1. Lecture Notes on Advanced Statistical Physics](#1-lecture-notes-on-advanced-statistical-physics)
  - [1.1 Recap: Equilibrium statistical mechanics](#11-recap-equilibrium-statistical-mechanics)
  - [Key concepts](#key-concepts)
  - [1.2 States and ensembles](#12-states-and-ensembles)
  - [1.3 Entropy and information](#13-entropy-and-information)
  - [1.4 Equilibrium ensembles](#14-equilibrium-ensembles)
  - [1.5 Averages, fluctuations and large deviations](#15-averages-fluctuations-and-large-deviations)
  - [1.6 Gartner-Ellis theorem](#16-gartner-ellis-theorem)
    - [1.6.1 Generic properties of $I(a)$](#161-generic-properties-of-ia)
    - [1.6.2 Interpreting GE in terms of a cumulating generating function](#162-interpreting-ge-in-terms-of-a-cumulating-generating-function)
      - [Example Gaussian:](#example-gaussian)
      - [Examples Gaussian Mixture](#examples-gaussian-mixture)
    - [1.6.3 Informal proof of GE-theorem](#163-informal-proof-of-ge-theorem)

---


## 1.1 Recap: Equilibrium statistical mechanics

Liquids, Solids, traffic jams, markets, societies, single molecule

Large Systems
:   Particle number $N =  N_a \in O(e23)$
    ==Analytically impossible== (⚡️ sometimes integratable systems)
    Simulations: $N \in O(e12)$ molecule dynamics

## Key concepts

Ensembles
:   ==Probability distribution== (ensembles)
    variances, ==fluctuation==
    determined by macroscopic variables $t, V, N$

Equilibrium
:   ingnores time variations
    slow changes (==thermodynamiccly reversable==)
    $\rightarrow$ Later focus on dynamics, non equilibrium

Statistical mechanics
:   derive macroscopic behavior from microscopic details, ==fluctuation== is important (eg. Variance of Energy)

Thermodynamics
:   $N → ∞$, fluctuations ignored, looking at averages
    general relations between macroscopic variables

## 1.2 States and ensembles

Microscopic state in class. mech.
:   Fully defined by position and momentum, ($(\vec{r}, \vec{p}) \in \R^{6N}$)

Ensembles
:   ==Probability distribution==
determined by macroscopic variables $t, V, N$
called ==macrostates==
$P(\vec{r},\vec{p}) = \frac{1}{N!} \frac{\rho(\vec{r},\vec{p})}{(2\pi\hbar)^{3N}}$

Average of generic observable $A(\vec{r},\vec{p})$
:   $\braket A  = \int d\vec{r},d\vec{p} P(\vec{r},\vec{p}) A(\vec{r},\vec{p}) \\
= \frac{1}{N!} \int \frac{d\vec{r},d\vec{p}}{2\pi\hbar} \rho(\vec{r},\vec{p}) A(\vec{r},\vec{p}) = Tr(A\rho)$

> To describe a QM Ensemble, use density operatir/matrix

In one normalized quantum state $ \ket i $
: $\ket i $ gives $\braket A  = \braket{i|A|i}$

QM average of observable distribution $p_i$ across $I$ states
:   $\braket A  = \sum_{i=1}^{I}p_i \braket{i|A|i} $
    $= \sum_{\mu}\sum_{i=1}^{I}p_i \braket{i|\mu}  \braket{\mu|A|i} $
    with $\mu$ an orthonormalbasis
    $= \sum_{\mu}\sum_{i=1}^{I}p_i \braket{\mu|A|i} \braket{i|\mu} $
    $= Tr(A\rho)$ with $\rho = \sum_{i=1}^{I}p_i \ket i \bra i $
    $\rho$ is the ==density matrix==

Properies of the density matrix $\rho$
:   Hermitian, non-negative, $Tr(\rho) = 1$
    Non-negativity: $\braket{\psi|\rho|ψ} = \sum_{i=1}^I p_i \braket{ψ|i} \braket{i|ψ} \geq 0$
    $Tr(\rho) = \sum_{μ} \braket{μ|ρ|μ} = ∑_{i=1}^{I} p_i \braket{μ|i}braket{i|μ} = ∑_{i=1}^I p_i=1$

What are Eigenvalues $\rho_n$ of $ρ$ and Eigenstates $\ket n$
:   $\rho_n \geq 0 $, $∑_n = 1$
    → $ρ_n$ playes ==the role of classical probabilities==
    → Pure state $\ket i$, $\rho_i = 1 $ all others 0

> discrete states in fixed volume because countable, thus $∑$ ist a good thing.

## 1.3 Entropy and information

> State $n$ with $\rho_n$ has ==information content== $-k_b \ln ρ_n$
~ *Shannon, source coding theorem*

Entropy $S$
:   Is the average information content.
$S = -∑_n k_b ρ_n \ln ρ_n = -k_b Tr(\rho \ln \rho)$ (QM)
$S = -k_b \frac{1}{N!} ∫\frac{d\vec{r} d\vec{p} }{(2\pi\hbar)^{3N}} ρ(\vec{r},\vec{p}) \ln ρ((\vec{r},\vec{p})$

Properties of Entropie $S$
:   $S ≥ 0$
    $S ≤ k_b \ln(\text{dim of state space})$ eg. for Spin: $S ≤ k_b N \ln 2$

> For proving bounds can use KL-divergence

KL-divergence
:   A measure of distance between probabillity distributions
    $KL(\rho'||ρ) = ∑_n ρ_n'\ln(ρ_n'/ρ_n) ≥ 0$

## 1.4 Equilibrium ensembles

The fundamental postulate of stat. mech.
:   In a closed system specified by energy $E$, particle number $N$ and volume $V$, ==all states are equally likely==

Microcanonical ensemble
:   fixed volume and particle number $V$, $N$ and energy $ E-Δ ≤ H ≤ E $ , then
    $p ∝ Θ(E-H)Θ(H-(E-Δ))$ ($Θ$ is Heaviside function)
    Maximizes entropy given constraints, Entropy is $S(E,N,V)$

> Mostly want to look at ==open== systems, which can exchange energy with "reservoir". If reseroir is large enough so that energy $\approx$ const.

Canonical Ensemble / Boltzmann distribution
:   energy $\approx$ const.
    partition function $Z = Tr{e^{-\beta H}}$
    distribution: $\rho=\frac{1}{Z}e^{-\beta H}$, where $\beta=1/k_BT$
    free energy: $F:= -k_BT\ln Z$
    $\rho$ also minimises free energy:
    $F[\rho']=\underbrace{Tr(H\rho')} - T\underbrace{(-k_B Tr(\rho')) \ln \rho'} \\ \qquad \qquad   \langle H \rangle ' \qquad \qquad \qquad S'$

Boltzmann const. $k_B$
:   $k_B = 1.380649 jK^{-1}$
    useful for estimates: $k_BT\approx0.025 eV$ at room temp.
    For theory often set to one

Grand canonical ensamble
:   If system can exchange ==energy== and ==particles== with reservoir
    density: $\rho = \frac{1}{Y}e^{-\beta (H-\mu N)}$, with $\mu$ chemical Potential and $Y=Tr(e^{-\beta(H_\mu N)})$
    Grandcanonical free energy: K = -k_BT\lnn (Y)

> Similar ensembles for volume fluctuations etc.

## 1.5 Averages, fluctuations and large deviations

Generic observable (extensive)
:   e.g. $A=M=\sum_{i=1}^N \sigma_i^Z$
    ==Average== can be obtained as
    $\braket A = Tr(A\rho) or \braket{A} = \frac{-\partial F}{\partial h}$ with $F=-k_BT\ln{Z}$ and $Z=Tr(e^{-\beta(H-hA)})$

Proof
:   Assume $[H,A] = 0$  Then $F=-k_BT\ln{\sum_n e^{-\beta(E_n-kA_n)}}$
    $\frac{\partial F}{\partial k} = -k_BT\frac{\sum_n \beta A_n e^{}-\beta(E_n-kA_n) }{\sum_n e^{-\beta(E_n-kA_n)}} = -\braket{A} $

$\square$

Fluctuation-dissipation theorem (FDT) static derivation
:   Take the second derivative of $F$:
    $\frac{\partial^2F}{\partial k^2} = \frac{\partial}{\partial k} \frac{\sum_n A_n e^{}-\beta(E_n-kA_n) }{\sum_n e^{-\beta(E_n-kA_n)}} = -\braket{A} ) $
    $= \frac{-\sum_n \beta A_n^2 e^{-\beta(...)}}{\sum_n e^{-\beta(...)}}-\frac{\sum_n \beta A_n e^{-\beta()}}{(\sum_ne^{-\beta(...)})}(-\sum_n A_n e^{-\beta(...)})$
    $=-\beta \braket{A^2}+\beta \braket{A}^2 = -\beta(\braket{A^2}- \braket{A}^2)=-\beta\braket{(\Delta A)^2}$,
    with $\Delta A$ as $A$s variancce

Static FDT
:   Relation of fluctuations and reactions on external currents.
    When $A$ is extensive, expect $F(k)$ also to be extensive ie.
    $F(k) = N f(k)$ for large $N$
    Then $\braket{A}=\frac{-\partial F}{\partial k}$ is $~N$, $\braket{(\Delta A)^2}\propto \frac{\partial^2 F}{\partial k^2}~N$
    So $\frac{\braket{(\Delta A)^2}^{1/2}}{\braket{A}}=$ relative standard deviation ~$\frac{N^{1/2}}{N}=N^{-1/2}\rightarrow 0$

> Typicall fluctuations $\Delta A ~ \sqrt{N}$ are Gaussia.

Gaussian distribution
: $P(a) = \frac{1}{\sqrt{2\pi \sigma^2/N }} \exp{(\frac{-N}{2\sigma}(a-a*)^2})$

## 1.6 Gartner-Ellis theorem

Rate function for $h=0$ is computes as
: $ I(a) = \sup_h \beta[ha + f(h)-f(0)] $, where $f(h)-f(0) = \delta_f(h)$

Extremum condition
:   $a= - \frac{\partial f}{\partial h} = \frac{\braket{A}_h}{N}$

### 1.6.1 Generic properties of $I(a)$

1. $I(a) \geq 0$
1. $I(a^*) = 0$, where $a^*$ mean value of $H, h=0$
1. $I(a)$: convex function
1. Quadratic expansion around $a^*$ gives Gaussian approx 
   $I(a)\equiv I(a^*)+ (a-a^*) I'(a^*) + \frac{(a-a^*)^2}{2}I''(a^*)$
   $=\frac{(a-a^*)^2}{2}I''(a^*)$

### 1.6.2 Interpreting GE in terms of a cumulating generating function

$$ \delta_f(h) = f(h)-f(0) = \frac{F(h)-F(0)}{N}  $$
$$ = -\frac{1}{N\beta}[\ln(Z(h))- \ln(Z(0))] $$
$$ = -\frac{1}{N\beta} \ln{\frac{Z(h)}{Z(0)}} $$
$$ = -\frac{1}{N\beta}\ln{\frac{Tr(e^{-\beta H} e^{\beta h A} )}{Tr(e^{-\beta H})}}$$
$$ = \frac{1}{N\beta} \ln{\braket{e^{\beta h A}}}_{|h=0} $$

$$ \frac{\partial^n}{\partial(\beta h)^n} \ln{\braket{e^{\beta h A}}}_{|h=0} = K_n $$, the $n$-th cumulant

$$ K_1 = \frac{\partial}{\partial(\beta h)} \ln{\braket{e^{\beta hA}}} = \frac{1}{\braket{e^{\beta hA}}} \braket{A e^{\beta h A}}_{|h=0} = \braket{A}$$

As Homework:
$$ K_2 = \braket{A^2} - \braket{A}^2 = Var{A} $$

#### Example Gaussian:

$A$ is gaussian distributed $P_G(A) = N(A|N_{a^*}, N_{\sigma^2})$
  
$$ \delta_f(h) = -\frac{1}{\beta N} \ln\braket{e^{\beta h A}} $$
$$ = -\frac{1}{\beta N} \int dA P_G(A) e^{\beta h A} $$
Homework:
$$ = -\frac{1}{\beta N} \ln e^{\beta h \braket{A} + (\beta h)^2/2 \braket{\Delta A^2}} $$
$$ = -\frac{1}{\beta N} \ln e^{\beta h N_{a^*} + (\beta h)^2/2 N_{\sigma^2}} $$
$$ \delta_f(h) = -ha^* - \frac{\beta h^2}{2} \sigma^2$$

Plugging this into GE-theorem:

$$ I(a) = \sup_h \beta[ha-ha^*-\frac{\beta h^2}{2}\sigma^2] $$

$$ \frac{\partial [...]}{\partial h} \overset{!}{=} 0 \Rightarrow h=\frac{(a-a^*)}{\beta \sigma^2}$$

Plugging result for $h$ into $I(a)$:

$$ I(a) = \beta \frac{(a-a^*)}{\beta \sigma^2} - \frac{\beta^2 (a-a^*)^2}{2\beta^2 \sigma^4}\sigma^2 $$

$$ = \frac{(a-a^*)^2}{2\sigma ^2} = I_G(a) $$

#### Examples Gaussian Mixture

$$ P(A) = \frac{1}{2} N(A|N_{a^*}, N_{\sigma^2}) + \frac{1}{2} N(A|-N_{a^*}, N_{\sigma^2}) $$
$$ P(a) =  \frac{1}{2} N(a|a^*, \sigma^2/N) + \frac{1}{2} N(a|-a^*, \sigma^2/N)) $$
$$ \tilde{I}(a) = -\lim_{N\rightarrow\infty} \frac{1}{N} \log P(a)$$
$$ \tilde{I}(a) = \min{\frac{(a-a^*)^2}{2\sigma^2},\frac{(a+a^*)^2}{2\sigma^2}}$$

Using GE:

$$ \delta_f(h) = -\frac{1}{\beta N} \ln\braket{e^{\beta h A}} $$
$$ =  -\frac{1}{\beta N} \ln \frac{1}{2}e^{\beta h N a^* + \frac{1}{2} \beta^2 h^2 N \sigma^2} + \frac{1}{2}e^{-\beta h N a^* + \frac{1}{2} \beta^2 h^2 N \sigma^2 } $$

Two cases:

$h>0$:
$$ \delta_f(h) =  -\frac{1}{\beta N} \ln e^{\frac{1}{2}\beta h N_a^*+ \frac{1}{2} \beta^2 h^2 N \sigma^2 + \text{exp. small terms}} $$
$$ = -ha^* - \frac{1}{2}\beta h^2 \sigma^2 $$

$h<0$:
$$\delta_f(h) = ha^* - \frac{1}{2}\beta h^2 \sigma^2$$

All in all:

$$ \delta_f(h) = -|h| a^* - \frac{1}{2} \beta h^2 \sigma^2$$

Because of $|h|$ we cant take a partial derivative but have to evaluate the supremum:

$$ I(a) = \sup_h \beta[ha - |h| a^* - \frac{1}{2} \beta h^2 \sigma^2]$$
cases:

- $a=0$:  $\sup_h[-\frac{1}{2}\beta h^2\sigma^2] = 0$ 
- $a=a^*$: $ \sup_h \beta[(h - |h|) a^* - \frac{1}{2} \beta h^2 \sigma^2] = 0$
- $a=-a^*$: $I(-a^*)=0$

As $\tilde{I}(a)$ is not convex this GE-yields its convex hull. Indicates phase transitions.

### 1.6.3 Informal proof of GE-theorem

(using saddle point approx)

$$ P(a) = \braket{\delta(\frac{A}{N}-a)} = N\braket{\delta(A-Na)}$$
expanding dirac $\delta$ with fourier
$$ = N\braket{\int \frac{d\lambda}{2\pi} e^{i\lambda(A-Na)}}$$
$$ = N \int  \frac{d\lambda}{2\pi} e^{-i\lambda Na} \braket{e^{i\lambda A}}$$

Given that:

$$ \braket{e^{i\lambda A}} = \frac{Tr(e^{-\beta H} e^{i\lambda A})}{Tr(e^{-\beta H})} = \frac{Z(i\lambda/\beta)}{Z(0)}$$

It follows:
$$ P(a) =  N \int \frac{d\lambda}{2\pi} e^{N[-i\lambda a - \beta f(i \lambda/\beta)+ \beta f(0)]}$$
Using saddle point approx. and only exponenial order:
$$ P(a) =   e^{N \text{extr}_h[-i\lambda a - \beta f(i \lambda/\beta)+ \beta f(0)]}$$
Set $h = i\lambda/\beta$
$$I(a) =\text{extr}_h[-\beta h a - \beta \delta_f(h)] = \text{extr}_h[\beta h a + \beta \delta_f(h)]$$
