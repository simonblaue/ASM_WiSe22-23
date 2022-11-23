---
title: "Lecture Notes on Advanced Statistical Physics 02 "
header-includes:
- \usepackage{braket}

output:
  pdf_document: 
       latex_engine: xelatex

---

- [2. Deterministic Dynamics](#2-deterministic-dynamics)
- [2.1 Quantum dynamics](#21-quantum-dynamics)
- [2.2 Classical dynamics](#22-classical-dynamics)
  - [Liouville equation](#liouville-equation)
  - [Correlation function](#correlation-function)

## 2. Deterministic Dynamics

- closed systems (no heat bath)

## 2.1 Quantum dynamics

Take Schrödinger and transform it into density matrix $\rho = \sum_n \rho_n \ket n \bra n$ has time evoultution (using $\frac{d}{dt}\ket{n} = (-i/\hbar) H \ket{n}$)

$$
\dot\rho = \sum_n \rho_n (\ket{\dot n}\bra{n} + \ket{n}\bra{\dot n}) \\
= \sum_n \rho_n (-\frac{i}{\hbar}H \ket{n}\bra{n}(\frac{i}{\hbar}H \ket{n}\bra{n})
$$

So $\dot\rho = -\frac{i}{\hbar} (H\rho - \rho H) =-\frac{i}{\hbar}[H,\rho]$
the ==von Neumann equation==

> Remember : $[H,\rho]$ is the commutator.

In equilibrium, $\rho$ is function of $H$ (microcanonical/Canonical)

$\Rightarrow$ $\dot\rho = 0$

Dynamics of averages
:   $\frac{d}{dt}\braket{A}= \frac{d}{dt}Tr(A\rho)= -\frac{i}{\hbar}Tr(A(H\rho-\rho H))$
    $= -\frac{i}{\hbar} Tr(AH\rho-HA\rho)=-\frac{i}{\hbar}Tr([A,H]\rho)$
    $ =  -\frac{i}{\hbar} \braket{[A,H]}$

> If you have a observable that commutes with $H$ the observed average of this quantatie is zero. In other words ==Conserved quantaties commute with $H$==

Formal solution of von Neumann equation
:   $\rho(t) = e^{-\frac{i}{\hbar}Ht}\rho(0)e^{\frac{i}{\hbar}Ht}$
    (Check: $\frac{d}{dt}\rho(t) = -\frac{i}{\hbar}H \underbrace{e^{...}\rho(0)e^{...}}_{\rho(t)}+\frac{i}{\hbar} \underbrace{e^{...}\rho(0)e^{...}}_{\rho(t)} H =  -\frac{i}{\hbar} (H\rho-\rho H)$)

$$ \braket{A(t)} = Tr(A\rho(t)) = Tr(Ae^{-\frac{i}{\hbar}Ht}\rho(0)e^{\frac{i}{\hbar}Ht}) \\
= Tr(\underbrace{e^{-\frac{i}{\hbar}Ht}Ae^{\frac{i}{\hbar}Ht}}_{A(t)}\rho(0))\\
$$

In the Heisenberg picture.

Representations of the evolution
:   1. $\rho(t)$ with fixed $A$
    1. $A(t)$ with fixed $\rho(0)$

Time evolution of $A(t)$
:   $\frac{d}{dt}A(t) = \frac{i}{\hbar}(HA(t)-A(t)H)= \frac{i}{\hbar}[H,A(t)]$

We can write formally:

$$
\frac{d}{dt}A(t) = L A(t) \\
\frac{d}{dt}\rho(t) = L^\dagger\rho(t)
$$

with $L(\cdot)= \frac{i}{\hbar}[H,\cdot]$ and $L^\dagger(\cdot)= -\frac{i}{\hbar}[H,\cdot]$

$L$ and $L^\dagger$ are adjoints of each other with respect to "product" defined as $Tr(AB)$, i.e. $Tr((L^\dagger A) B) = Tr(A (LB))$

> analog of transposing a matrix

Adjointness follows, because ==averages must be independent of the picture used==.

$$
\frac{d}{dt}\braket{A(t)} = Tr(A\dot\rho(t)) = Tr(A(L^\dagger \rho(t)))\\
\overset{!}{=} Tr(\dot A(t)\rho(0)) = Tr((LA(t)) \rho(0))
$$
At $t=0$, where $A(0)=A$, we get
$$
Tr(A (L^\dagger \rho(0))) = Tr((LA)\rho(0))
$$

## 2.2 Classical dynamics

$N$ particles, coordinates $r_{i\alpha}$, $i$ for number of particels, $\alpha$ for dimensions. Momenta $p_{i\alpha}$

Observable $A(\underline{r},\underline{p})$. What is the meaning of and how is $A(\underline{r}, \underline{p}, t)$ defined? We want:

$$
\braket{A(t)}= Tr(A(t)\rho(0)) = \frac{1}{N!}\int\frac{d\underline{r}d\underline{p}}{(2\pi\hbar)^Nd} A(\underline{r}, \underline{p}, t) \rho(\underline{r}, \underline{p}, 0)
$$

$A(\underline{r}, \underline{p}, t) = $ value of $A$ at time $t$ when starting from $\underline{r}, \underline{p}$

Time evolution of $A(\underline{r}, \underline{p}, t)$
:   $A(\underline{r}, \underline{p}, t+dt) = A(\underline{r}+d\underline{r}, \underline{p}+d\underline{p}, t)$
    $= A(\underline{r}, \underline{p}, t)+ \sum_{i\alpha}(\frac{\partial A}{\partial r_{i\alpha}} dr_{i\alpha}+ \frac{\partial A}{\partial p_{i\alpha}} dp_{i\alpha}) + ...$
    $\frac{A(\underline{r}, \underline{p}, t+dt)- A(\underline{r}, \underline{p},t)}{dt} = \sum_{i\alpha}(\frac{\partial A}{\partial r_{i\alpha}} \frac{dr_{i\alpha}}{dt}+ \frac{\partial A}{\partial p_{i\alpha}} \frac{dp_{i\alpha})}{dt}$
    $\Rightarrow \frac{\partial A}{\partial t} = \sum_{i\alpha}(\frac{\partial A}{\partial r_{i\alpha}} \frac{\partial H}{\partial p_{i\alpha}}+ \frac{\partial A}{\partial p_{i\alpha}} \frac{\partial H}{\partial r_{i\alpha})} ) =  \{A,H\} $

> Now oberve this in other picture, fixed observables -> How does distribution evolv?

As in Quantum Mechanics, we can represent time evolution alternatively in terms of $\rho(\underline{r}, \underline{p}, t)$ From the adjointness argument we get:

### Liouville equation

$\frac{\partial \rho}{\partial t} = - \sum_{i\alpha} [\frac{\partial }{\partial r_{i\alpha}} (\frac{\partial H}{\partial p_{i\alpha}} \rho) - \frac{\partial }{\partial p_{i\alpha}} (\frac{\partial H}{\partial r_{i\alpha})} \rho)]$
$= - \underline{\nabla} \cdot(\underline{v}\rho)$, which is the the continuity equation in phase space

$\underline{\nabla} = (\frac{\partial}{\partial r_{11}}, ... , \frac{\partial}{\partial p_{11}}, ...), \underline{v} = (\frac{\partial H}{\partial p_{11}}, ... , \frac{\partial H}{\partial r_{11}}, ...)$

For Hamiltonian dynamics, phase space flow is incompressible $(\underline{\nabla} \cdot \underline{v}) = 0$ so equivalently:

$$ \frac{\partial \rho}{\partial t} = -\underline{v}\cdot \underline{\nabla} \rho = - \sum_{i\alpha} [ \frac{\partial H}{\partial p_{i\alpha}} \frac{\partial \rho}{\partial r_{i\alpha}} -  \frac{\partial H}{\partial r_{i\alpha})} \frac{\partial \rho}{\partial p_{i\alpha}} ] = \{H,\rho\} = -\{\rho, H\}$$

We can write:

$\frac{\partial A}{\partial t} = LA, L(\cdot) = \{\cdot,H\}$

$\frac{\partial \rho}{\partial t} = L^\dagger A, L(\cdot) = \{H,\cdot\}$

These are ==adjoints== with respect to $\int d\underline{r}d\underline{p}A(\underline{r},\underline{p}) B(\underline{r},\underline{p})$

Derivation of Liouville equation
:   we will omit the prefactor $\frac{1}{N!}\frac{1}{(2\pi\hbar)^{3N}}$
    Equallity of $\braket{A(t)}$ in the two pictures gives
    $\int d\underline{r}d\underline{p} \rho(\underline{r},\underline{p}, 0) A(\underline{r},\underline{p}, t) = \int d\underline{r}d\underline{p} \rho(\underline{r},\underline{p}, t) A(\underline{r},\underline{p})$

Take $\frac{d}{dt}$ at $t=0$:

$$
\int  d\underline{r}d\underline{p} \frac{\partial \rho}{\partial t} A = \int d\underline{r}d\underline{p} \rho(0) \frac{\partial A}{\partial t}\\
= \int d\underline{r}d\underline{p} \rho(0) \{A,H\} \\
= \int d\underline{r}d\underline{p} \rho(0) \sum_{i\alpha}(\frac{\partial A}{\partial r_{i\alpha}} \frac{\partial H}{\partial p_{i\alpha}}+ \frac{\partial A}{\partial p_{i\alpha}} \frac{\partial H}{\partial r_{i\alpha})} )
$$

Integrate by parts:

$= \int d\underline{r}d\underline{p} \sum_{i\alpha} A(\frac{\partial}{\partial r_{i\alpha}}(\rho(0)\frac{\partial H}{\partial p_{i\alpha}}) +\frac{\partial}{\partial p_{i\alpha}}(\rho(0)\frac{\partial H}{\partial r_{i\alpha}})$

Must be true for ==all== $A$, so

$$
\frac{\partial \rho}{\partial t} = -\sum_{i\alpha} \frac{\partial}{\partial r_{i\alpha}}(\rho \frac{\partial H}{\partial p_{i\alpha}}) -\frac{\partial}{\partial p_{i\alpha}}(\rho\frac{\partial H}{\partial r_{i\alpha}}
$$

So far: $\frac{\partial A}{\partial t} = LA$, $LA=-\{H,A\}$ or $= \frac{i}{\hbar}[H,A]$

Formal solution: $A(t)= e^{Lt}A$

### Correlation function

Define ==correlation function==
:   $C_{AB}(t,t') = Tr A(t) B(t') \rho(0)$, two-time function in general
    Symmetries:
    $C_{AB}(t,t') = C_{BA}(t,t')$ in classical case
    $C^*_{AB}(t,t') = C_{BA}(t,t')$ in quantum case take conjugent.

Focus now on equilibrium relations $\rho(0) = \rho^{\text{eq}}$. Then we can write the correlation function as a scalar product:

$ C_{AB}(t,t')= (A(t),B(t'))$ with  $(A,B)= Tr(AB\rho^{\text{eq}})$

$C_{AB}(t,t') = (e^{\mathscr{L}t} A,e^{\mathscr{L}t'} B)$

Can show that $\mathscr{L}$ is anti-selfadjoint w.r.t. $(A,B)$: $(A,\mathscr{L}B) = - (\mathscr{L}A, B)$

Proof QM version
:   exploit that $[H,\rho^{\text{eq}}]= 0$
    $(A,\mathscr{L}B) = \frac{i}{\hbar}(A,[H,B])= \frac{i}{\hbar}Tr(A(HB-BH)\rho^{\text{eq}})$
    $= \frac{i}{\hbar}Tr(AHB\rho^{\text{eq}}-AB\rho^{\text{eq}}H)$
    $= \frac{i}{\hbar}Tr(AHB-HAB)\rho^{\text{eq}} = -\frac{i}{\hbar}(HA-AH,B) = -(\mathscr{L}A,B)$

Consequence
:   $(A,\mathscr{L}^n B) = (-1)^n(\mathscr{L}A,B)$
    $(A,e^{\mathscr{L}t'}B) = \sum_{n=0}^{\infin} \frac{1}{N!}(t')^n(A,\mathscr{L}^nB) = (e^{-\mathscr{L}t'}A,B)$
    So $C_{AB}(t,t') = (e^{\mathscr{L}t}A,e^{\mathscr{L}t'}B) = (e^{-\mathscr{L}t'}e^{\mathscr{L}t}A,B) = (e^{-\mathscr{L}(t-t')}A,B) = C_{AB}(t-t')$ as expected in equilibrium

For auto-correlation $(A=B): C_{AA}^*(t,t')=C_{AA}(t,t')$ so in equilibirum
$C_{AA}^*(\underbrace{t-t'}_{\Delta t}) = C_{AA}(\underbrace{t'-t}_{-\Delta t})$

Classical case: $C_{AA}(\Delta t) = C_{AA}(-\Delta t)$ symmetric in time
Can solve that $|C_{AA}(\Delta t)| \leq C_{AA}(0)$
