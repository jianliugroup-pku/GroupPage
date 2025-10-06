# A general suite for benchmarking nonadiabatic dynamics

[TOC]

In this tutorial we provide a general suite for benchmarking nonadiabatic dynamics described in ref [^1]  together with [data](files/data.zip) (this link is the collection) of all figures. Among all the dynamics tested, the robustness of nonadiabatic field (NaF)[^a][^b] is highlighted.

For a full publication list of the generalized constraint coordinate-momentum phase space (CPS) and NaF please see [this page](publication_list.html).

For an illustration of the nonadiabatic force on the surface hopping dynamics please see [this page](fsnaf.html) (reported in ref[^a]).

## 1. LVCM (Linear vibronic coupling models)

[Figures_and_data](LVCM.html)

The Hamiltonian of a LVCM (linear vibronic coupling model) reads
$$
\hat H = {\hat H_0} + {\hat H_C},
$$
where
$$
{\hat H_0} = \sum\limits_{k = 1}^{{N_{{\rm{nuc}}}}} {\frac{{{\omega _k}}}{2}\left( {\hat {\bar{P}}_k^2 + \hat{\bar{R}}_k^2} \right)},
$$

$$
\hat{H}_C = \sum^F_{n=1}\left(E_n+\sum^{N_\textrm{nuc}}_{k=1}\kappa^{(n)}_k \hat{\bar{R}}_k\right)\left\vert n\right\rangle\left\langle n\right\vert+\sum^F_{n\neq m}\left(\sum^{N_\textrm{nuc}}_{k=1}\lambda_k^{(nm)}\hat{\bar{R}}_k\right)\left\vert n\right\rangle\left\langle m\right\vert.
$$

Note that the relation between the dimensionless weighted normal-mode coordinate/momentum, $\{\bar{R}_k,\bar{P}_k\}$, and the corresponding canonical (mass-weighted) coordinate/momentum, $\{R_k,P_k\}$, in the diabatic representation reads
$$
\bar{R}_k = \sqrt{\omega_k} R_k,\quad \bar{P}_k = P_k/\sqrt{\omega_k}
$$
The parameters of LVCMs in the suite are listed as follows, except the parameters of the 7-state 39-mode LVCM of Thymine in ref [^2] which are kindly provided by Martha Yaghoubi Jouybari and Fabrizio Santoro.

**Table 1. Parameters of 2-state 3-mode LVCM of Pyrazine (unit: eV) **[^3]

| Parameters                                                   | Value                        |
| ------------------------------------------------------------ | ---------------------------- |
| $E_1,E_2$                                                    | 3.94, 4.84                   |
| $\kappa_1^{(1)}, \kappa_2^{(1)}, \kappa_1^{(2)}, \kappa_2^{(2)}$ | 0.037, -0.105, -0.254, 0.149 |
| $\lambda_{3}^{(12)}$                                         | 0.262                        |
| $\omega_1, \omega_2, \omega_3$                               | 0.126, 0.074, 0.118          |

**Table 2. Parameters of 24-mode LVCM of Pyrazine (unit: eV)**[^4]

| Parameters           | Value           |
| -------------------- | --------------- |
| $E_1, E_2$           | -0.4617, 0.4617 |
| $\lambda_{1}^{(12)}$ | 0.1825          |

| Mode | $\omega$ | $\kappa^{(1)}$ | $\kappa^{(2)}$ |
| ---- | -------- | -------------- | -------------- |
| 1    | 0.0936   |                |                |
| 2    | 0.074    | -0.0964        | 0.1194         |
| 3    | 0.1273   | 0.0470         | 0.2012         |
| 4    | 0.1568   | 0.1594         | 0.0484         |
| 5    | 0.1347   | 0.0308         | -0.0308        |
| 6    | 0.3431   | 0.0782         | -0.0782        |
| 7    | 0.1157   | 0.0261         | -0.0261        |
| 8    | 0.3242   | 0.0717         | -0.0717        |
| 9    | 0.3621   | 0.0780         | -0.0780        |
| 10   | 0.2673   | 0.0560         | -0.0560        |
| 11   | 0.3052   | 0.0625         | -0.0625        |
| 12   | 0.0968   | 0.0188         | -0.0188        |
| 13   | 0.0589   | 0.0112         | -0.0112        |
| 14   | 0.0400   | 0.0069         | -0.0069        |
| 15   | 0.1726   | 0.0265         | -0.0265        |
| 16   | 0.2863   | 0.0433         | -0.0433        |
| 17   | 0.2484   | 0.0361         | -0.0361        |
| 18   | 0.1536   | 0.0210         | -0.0210        |
| 19   | 0.2105   | 0.0281         | -0.0281        |
| 20   | 0.0778   | 0.0102         | -0.0108        |
| 21   | 0.2294   | 0.0284         | -0.0284        |
| 22   | 0.1915   | 0.0196         | -0.0196        |
| 23   | 0.4000   | 0.0306         | -0.0306        |
| 24   | 0.3810   | 0.0269         | -0.0269        |

In our simulation of the pyrazine models, the second diabatic (electronic) state is initially occupied, and the nuclear initial condition is set to be the Wigner distribution of the vibrational ground state,
$$
\rho_W(\mathbf{\bar{R},\bar{P}})\propto \exp\left[-\sum^N_{k=1}\left(\bar{R}_k^2+\bar{P}_k^2\right)\right].
$$
The numerically exact results in our data are produced by MCTDH[^5].

**Table 3. Parameters of 3-state 2-mode LVCM of $\mathrm{Cr(CO)_5}$ (unit: eV)**[^6]

| Parameters                                                   | Value                    |
| ------------------------------------------------------------ | ------------------------ |
| $E_1, E_2, E_3$                                              | 0.0424, 0.0424, 0.4344   |
| $\kappa_2^{(1)}, \kappa_2^{(2)}$                             | -0.0328, 0.0328          |
| $\lambda_{1}^{(12)}, \lambda_{1}^{(23)}, \lambda_{2}^{(13)}$ | 0.0328, -0.0978, -0.0978 |
| $\omega_1, \omega_2$                                         | 0.0129, 0.0129           |

In our test of the $\mathrm{Cr(CO)_5}$ model, the second diabatic (electronic) state is initially occupied, and the nuclear initial condition is set to be the Wigner distribution of the following Wigner distribution,
$$
\rho_W(\mathbf{\bar{R},\bar{P}})\propto\exp\left[-\sum^2_{k=1}\left(\frac{\left(\bar{R}_k-r\right)^2}{2\alpha_k^2}+2\alpha_k^2\bar{P}_k^2\right)\right],
$$
where $r_1 = 0, r_2=14.3514, \alpha_1 = 0.4501$, and $\alpha_2=0.4586$. The numerically exact results in our data are produced by MCTDH, taken from ref [^6].

In our test of the Thymine model, the $\pi\pi^*2$ state is initially occupied with each nuclear mode in its vibrational ground state, i.e., eq (5). The numerically exact results in our data are produced by ML-MCTDH, taken from ref [^2].

## 2. Tully models

[Figures_and_data](Tully.html)

The single avoided crossing (SAC) model[^7] reads
$$
\begin{aligned}
    &V_{11}(R)=-V_{22}(R)=A[1-\exp(-B|R|)]\text{sgn}(R),\nonumber\\
    &V_{12}(R)=V_{21}(R)=C\exp(-DR^2),
\end{aligned}
$$
where $A=0.01, B=1.6, C=0.005, D=1$ (in a.u.). The dual avoided crossing (DAC) model[^7] reads
$$
\begin{aligned}
    &V_{11}(R)=0, V_{22}(R)=-A\exp(-BR^2)+E_0, \nonumber\\
    &V_{12}(R)=V_{21}(R)=C\exp(-DR^2),
\end{aligned}
$$
with $A=0.01, B=1.6, C=0.005, D=1.0$ (in a.u.). The extended coupling region (ECR) model [^7] reads
$$
   \begin{aligned}
   &V_{11}(R)=-V_{22}(R)=E_0, \nonumber\\
   &V_{12}(R)=V_{21}(R)=C[\exp(BR)h(-R)+(2-\exp(-BR))h(R)],
   \end{aligned}
$$
where $E_0 = -0.0006, B=0.9, C=0.1$ (in a.u.), and $h(R)$ denotes the Heaviside step function. In our data, the nuclear DOF is with mass $m=2000$ (in a.u.). The electronic ground state is selected as the initial state, while the initial condition for the nuclear DOF follows the Wigner distribution
$$
\rho_W(R,P)\propto\exp\left(-\alpha(R-R_0)^2-(P-P_0)^2/\alpha\right)
$$
with width parameter $\alpha=1$. The coordinate center $R_0$ for SAC, DAC, and ECR models is set to -3.8, -10.0, and -13.0, respectively. The initial momentum $P_0$ is scanned and is used as the x-axis in our figures.

The exact results in our data are provided by DVR[^8].

## 3. 3-state photodissociation models

[Figures_and_data](Morse3.html)

The potential matrix elements of anharmonic 3-state photodissociation models[^9] reads
$$
\begin{aligned}
{V_{ii}}\left( R \right) &= {D_i}{\left[ {1 - {e^{ - {\beta _i}\left( {R - {R_i}} \right)}}} \right]^2} + {C_i},{\rm{  }}\;i = 1,2,3.\\
{V_{ij}}\left( R \right) &= {V_{ji}}\left( R \right) = {A_{ij}}{e^{ - {\alpha _{ij}}{{\left( {R - {R_{ij}}} \right)}^2}}}{\rm{,  }}\;\;\;i,j = 1,2,3;{\text{ and }}i \ne j.
\end{aligned}
$$
The parameters are listed as follows,

**Table 4. Parameters of the 3-state photodissociation models (unit: a.u.)**[^9]

| Parameters                              | Model 1           | Model 2           | Model 3           |
| --------------------------------------- | ----------------- | ----------------- | ----------------- |
| $D_{1}, D_{2}, D_{3}$                   | 0.003,0.004,0.003 | 0.020,0.010,0.003 | 0.020,0.020,0.003 |
| $\beta_{1}, \beta_{2}, \beta_{3}$       | 0.65,0.60,0.65    | 0.65,0.40,0.65    | 0.40,0.65,0.65    |
| $R_{1}, R_{2}, R_{3}$                   | 5.0,4.0,6.0       | 4.5,4.0,4.4       | 4.0,4.5,6.0       |
| $C_{1}, C_{2}, C_{3}$                   | 0.00,0.01,0.006   | 0.00,0.01,0.02    | 0.02,0.00,0.02    |
| $A_{12}, A_{23}, A_{31}$                | 0.002,0.002,0.0   | 0.005,0.0,0.005   | 0.005,0.0,0.005   |
| $R_{12}, R_{23}, R_{31}$                | 3.40,4.80,0.00    | 3.66,0.00,3.34    | 3.40,0.00,4.97    |
| $\alpha_{12}, \alpha_{23}, \alpha_{31}$ | 16.0,16.0,0.0     | 32.0,0.0,32.0     | 32.0,0.0,32.0     |
| $R_e$                                   | 2.9               | 3.3               | 2.1               |

The initial electronic configuration is set to the first diabatic state, and the nuclear DOF is sampled from the Wigner distribution of the ground state,
$$
\rho_W(R,P) \propto\exp\left(-m\omega(R-R_e)^2-P^2/m\omega\right),
$$
where $m=20000$ a.u. is the mass of the nuclear DOF, $\omega=0.005$ a.u. is the vibrational frequency of the ground state, and the center $R_e$ is listed in Table 4.

The exact results in our data are produced by DVR[^8].

## 4. 2-state photodissociation model

[Figures_and_data](Morse2.html)

The potential matrix elements of the anharmonic 2-state photodissociation model[^9][^1] reads
$$
\begin{aligned}
{V_{ii}}\left( R \right) &= {D_i}{\left[ {1 - {e^{ - {\beta _i}\left( {R - {R_i}} \right)}}} \right]^2} + {C_i},{\rm{  }}\;i = 1,2.\\
{V_{ij}}\left( R \right) &= {V_{ji}}\left( R \right) = {A_{ij}}{e^{ - {\alpha _{ij}}{{\left( {R - {R_{ij}}} \right)}^2}}}{\rm{,  }}\;\;\;i,j = 1,2;{\text{ and }}i \ne j.
\end{aligned}
$$
The parameters are listed as follows,

**Table 5. Parameters of the 2-state photodissociation model (unit: a.u.)**[^9][^1]

| Parameters             | Values      |
| ---------------------- | ----------- |
| $D_{1}, D_{2}$         | 0.020,0.003 |
| $\beta_{1}, \beta_{2}$ | 0.65,0.65   |
| $R_{1}, R_{2}$         | 4.5,4.4     |
| $C_{1}, C_{2}$         | 0.00,0.02   |
| $A_{12}$               | 0.005       |
| $R_{12}$               | 3.34        |
| $\alpha_{12}$          | 8.0         |
| $R_e$                  | 3.3         |

The initial condition and mass are identical to those of Model 2 of the 3-state photodissociation models.

The exact results in our data are produced by DVR[^8].

## 5. Spin-boson models

[Figures_and_data](spin_boson.html)

The Hamiltonian of a spin-boson model reads[^10]
$$
\begin{aligned}
    \hat{H}=\hat{H}_b+\hat{H}_s+\hat{H}_{b-s}=\sum_{j=1}^{N_b}\left(\frac{\hat P_j^2}{2}+\frac{1}{2}\omega_j^2\hat R_j^2\right)+
    \left[\begin{matrix}
    \varepsilon & \Delta \\
     \Delta & -\varepsilon 
    \end{matrix}\right]
    +\left[\begin{matrix}
    1 & 0\\
     0 & -1 \end{matrix}\right]\sum_{j=1}^{N_b}c_j\hat R_j  
\end{aligned}
$$
We use the discretization scheme for the Ohmic spectral density[^10] with the Kondo parameter $\alpha$ and the cut-off frequency $\omega_c$, which results in the following expressions[^11][^12] with $N_b$ discrete modes:
$$
    \omega_j=-\omega_c\ln\left(1-\frac{j}{1+N_b}\right),\quad  c_j=\omega_j\sqrt{\frac{\alpha\omega_c}{N_b+1}},\quad j=1,\cdots,N_b
$$
We employ four specific spin-boson models with $\varepsilon=\Delta=1$ at $\beta = 5$. The parameters[^30] of the bath DOFs are $\alpha\in\{0.1,0.4\}$ and $\omega_c\in\{1,2.5\}$ with $N_b=300$.

The system occupies the first diabatic electronic state, and the nuclear DOFs are sampled from their equilibrium Wigner density. The numerically exact results in our data[^30] are calculated by eHEOM[^13].

## 6. Site-exciton models

[Figures_and_data](site_exciton.html)

The Hamiltonian of site-exciton models reads[^31]
$$
 \hat{H}=\hat{H}_{s}+\hat{H}_{b}+\hat{H}_{s-b},  
$$
where $\hat{H}_{b}=\sum_{n, j} \frac{1}{2}\left(\hat{P}_{n j}^{2}+\omega_{n j}^{2} \hat{R}_{n j}^{2}\right)$ and $\hat{H}_{s-b}=\sum_{n, j} c_{n j} \hat{R}_{nj}|n\rangle\langle n|$. We employ the Debye spectral density
$$
    J(\omega)=\frac{2\lambda\omega_c\omega}{\omega_c^2+\omega^2}
$$
for each state, where $\lambda$ and $\omega_c$ denote the bath reorganization energy and the characteristic frequency, respectively. The corresponding discretization scheme[^14][^15][^16] is
$$
\begin{aligned}
&\omega_{nj}=\omega_c\tan\left[\frac{\pi}{2}\left(1-\frac{j}{N_b+1}\right)\right] \\
    &c_{nj}=\sqrt{\frac{2\lambda}{N_b+1}}\omega_{nj}, \\
    & n=1,\cdots, F; \quad j=1,\cdots, N_b
    \end{aligned}
$$
For the 7-state FMO model[^17], the system Hamiltonian reads
$$
\hat{H}_{s}=\left(\begin{array}{ccccccc}
12410 & -87.7 & 5.5 & -5.9 & 6.7 & -13.7 & -9.9 \\
-87.7 & 12530 & 30.8 & 8.2 & 0.7 & 11.8 & 4.3 \\
5.5 & 30.8 & 12210 & -53.5 & -2.2 & -9.6 & 6.0 \\
-5.9 & 8.2 & -53.5 & 12320 & -70.7 & -17.0 & -63.3 \\
6.7 & 0.7 & -2.2 & -70.7 & 12480 & 81.1 & -1.3 \\
-13.7 & 11.8 & -9.6 & -17.0 & 81.1 & 12630 & 39.7 \\
-9.9 & 4.3 & 6.0 & -63.3 & -1.3 & 39.7 & 12440
\end{array}\right) \mathrm{cm}^{-1},
$$
The bath reorganization energy is $\lambda=35 \mathrm{~cm}^{-1}$, the characteristic frequency is $\omega_{c}=106.14 \mathrm{~cm}^{-1}$, and $N_b = 50$ for each state. 

For 3-state singlet-fission model of pentacene[^18][^19], the system Hamiltonian reads
$$
\hat H_s=\left[\begin{matrix}
    200 & -50 & 0 \\
    -50 & 300 & -50 \\
    0 & -50 & 0
    \end{matrix}\right]\text{ meV}
$$
The bath reorganization energy is $\lambda=0.1$ eV, the characteristic frequency is $\omega_c=0.18$ eV, and $N_b = 200$ for each state.

In these two site-exciton models, the bath DOFs are sampled from their equilibrium Wigner distribution, and the first diabatic electronic state is initially occupied.

The numerically exact results in our data[^1] are produced by TD-DMRG[^20] in our FMO model, and by HEOM[^21] in our SF model.

## 7. Atom-in-cavity model

[Figures_and_data](aic.html)

The total Hamiltonian for the atom-in-cavity models[^22][^23] can be decomposed into three parts
$$
\hat H=\hat H_a + \hat H_p + \hat H_c
$$
where
$$
\hat H_a=\sum_{n=1}^F\varepsilon_n\ket{n}\bra{n}
$$

$$
\hat H_p=\sum_{j=1}^N\frac{1}{2}\left(\hat P_j^2+\omega_j^2\hat R_j^2\right)
$$

$$
\hat H_c=\sum_{n\ne m}^F\left(\sum_{j=1}^N\omega_j\lambda_j(r_0)\hat R_j\right)\mu_{nm}\ket{n}\bra{m}
$$

The quantity $\mu_{nm}=\mu_{mn}$ denotes the transition dipole moment between the $n$-th and $m$-th atomic energy levels. The coupling strength between the atom and the $j$-th optical field mode is given by
$$
\lambda_j(r_0)=\sqrt{2/\varepsilon_0L}\sin(j\pi r_0/L),
$$
where $L=236200$ a.u. represents the cavity length, $r_0=L/2$ corresponds to the cavity center, and $\varepsilon_0$ is the vacuum permittivity. The mode frequency is determined by $\omega_j=j\pi c/L$, with $c=137.036 \textrm{\ a.u.}$ to be the light speed in vacuum.  We employ 400 standing wave modes and model the hydrogen atom system as a three-level system. The energy levels and transition dipole moments (in atomic units) are: $\varepsilon_1=-0.6738$, $\varepsilon_2=-0.2798$, $\varepsilon_3=-0.1547$, $\mu_{12}=-1.034$, $\mu_{23}=-2.536$, and $\mu_{13}=0$.  A reduced two-level system is also considered, where only the lowest two energy levels of the hydrogen atom is retained.

Initially the highest atomic state is occupied, with each optical field mode in their optical vacuum state with Wigner distribution reading
$$
\rho_W(R_j,P_j)\propto\exp\left(-P_j^2/\omega_j-\omega_jR_j^2\right)
$$
The exact results are produced by truncated configuration interaction and are taken from refs[^22][^23].

## 8. Holstein model

[Figures_and_data](mobility.html)

The Hamiltonian of the 1-dim Holstein model[^24][^25] is given by

$$
\hat H=\hat H_s + \hat H_p + \hat H_c
$$
where the electronic system is described by a nearest-neighbor tight-binding model:  
$$
\hat H_s=\sum_n E_n\hat a_{n}^\dagger\hat a_{n} + V(\hat a_{n+1}^\dagger\hat a_{n}+\hat a_{n}^\dagger\hat a_{n+1})
$$
Here, $\hat a_n^\dagger$ and $\hat a_n$ denote the fermionic creation and annihilation operators at the $n$-th lattice site, while $E_n$ and $V$ represent the on-site energy and the coupling between neighboring sites, respectively. For a periodic chain with a finite number of lattice sites $L$, one typically sets $E_n=0$ and assumes $\hat a_{n+L}=\hat a_n$.  

The phonon Hamiltonian describes a bosonic bath environment:  
$$
\hat H_p=\sum_j\sum_n \hbar\omega_{j} (\hat b_{jn}^\dagger\hat b_{jn}+1/2)
$$
The term $\hat H_c$ encompasses all electron-phonon coupling effects:  
$$
\hat H_c=\sum_{n,j} c_j^{(1)}(\hat b_{jn}^\dagger+\hat b_{jn})\hat a_n^\dagger\hat a_n
$$
Considering the low carrier concentration regime, where the system contains only a single electron. In this case, the fermionic operators act within the single-particle Fock space, and the system Hamiltonian becomes isomorphic to a multistate Hamiltonian[^26]. Consequently, the fermionic operators can be expressed in terms of state representations:  
$$
\begin{aligned}

  &\hat a_n^\dagger=\ket{0_1,\cdots,1_n,\cdots,0_L}\bra{0_1,\cdots,0_n,\cdots,0_L}\triangleq \ket{n}\bra{\tilde 0},\\

  &\hat a_n=\ket{0_1,\cdots,0_n,\cdots,0_L}\bra{0_1,\cdots,1_n,\cdots,0_L}\triangleq \ket{\tilde 0}\bra{n}

\end{aligned}
$$
or equivalently,  
$$
\hat a_n^\dagger\hat a_m=\ket{n}\bra{m}
$$


We set $V = 0.083 \mathrm{\ eV}, $$d=R=7.2\mathrm{\ Å}$ and $L=21$[^27][^28], where the unmentioned parameters are used in the calculation of correlation function. Other parameters are listed as follows.

**Table 6. Parameters of the Holstein model in refs**[^27][^28]

| Mode | $\omega_j$ (meV) | $g_j=c_j^{(1)}/\hbar\omega_j$ |
| ---- | ---------------- | ----------------------------- |
| 1    | 10.0             | 0.96                          |
| 2    | 27.0             | 0.38                          |
| 3    | 78.0             | 0.25                          |
| 4    | 124.0            | 0.20                          |
| 5    | 149.0            | 0.15                          |
| 6    | 167.0            | 0.31                          |
| 7    | 169.0            | 0.13                          |
| 8    | 190.0            | 0.20                          |
| 9    | 198.0            | 0.31                          |

We guide readers to Section S5-S6 of Supporting Information of ref[^1] for details of initial condition and correlation function. The exact results in our data are taken from ref[^27], produced by TD-DMRG.

---

## 9. References

[^1]: B. Wu, B. Li, X. He, X. Cheng, J. Ren and J. Liu, "Nonadiabatic field: A conceptually novel approach for nonadiabatic quantum molecular dynamics ", J. Chem. Theory Comput. 21, 3775-3813 (2025).

[^2]: J. A. Green, M. Y. Jouybari, D. Aranda, R. Improta, and F. Santoro, “Nonadiabatic absorption spectra and ultrafast dynamics of DNA and RNA photoexcited nucleobases”, Molecules 26, 1743 (2021).
[^3]:R. Schneider and W. Domcke, “S1-S2 conical intersection and ultrafast S2→S1 internal-conversion in pyrazine”, Chem. Phys. Lett. 150, 235–242 (1988).
[^4]:S. Krempl, M. Winterstetter, H. Plöhn, and W. Domcke, “Path-integral treatment of multi-mode vibronic coupling”, J. Chem. Phys. 100, 926–937 (1994).
[^5]:G. A. Worth, M. H. Beck, A. Jackle and H.-D. Meyer, The MCTDH Package, Version 8.2, (2000). H.-D. Meyer Version 8.3 (2002), Version 8.4 (2007). O. Vendrell and H.-D. Meyer Version 8.5 (2013). Version 8.5 contains the ML-MCTDH algorithm. See http://mctdh.uni-hd.de. (accessed on November 1st, 2023) Used version: 8.5.14.
[^6]:G. A. Worth, G. Welch, and M. J. Paterson, “Wavepacket dynamics study of Cr(CO)5 after formation by photodissociation: Relaxation through an (E ⊕ A) ⊗ e Jahn–Teller conical intersection”, Mol. Phys. 104, 1095–1105 (2006).
[^7]:J. C. Tully, “Molecular-dynamics with electronic-transitions”, J. Chem. Phys. 93, 1061–1071 (1990).
[^8]: D. T. Colbert and W. H. Miller, "A novel discrete variable representation for quantum mechanical reactive scattering via the S-matrix Kohn method", J. Chem. Phys. 96, 1982-1991 (1992).
[^9]: E. A. Coronado, J. Xing, and W. H. Miller, “Ultrafast non-adiabatic dynamics of systems with multiple surface crossings: a test of the Meyer-Miller Hamiltonian with semiclassical initial value representation methods”, Chem. Phys. Lett. 349, 521–529 (2001).
[^10]:A. J. Leggett, S. Chakravarty, A. T. Dorsey, M. P. A. Fisher, A. Garg, and W. Zwerger, “Dynamics of the dissipative two-state system”, Rev. Mod. Phys. 59, 1–85 (1987).
[^11]:N. Makri, “The linear response approximation and its lowest order corrections: An influence functional approach”, J. Phys. Chem. B 103, 2823–2829 (1999).

[^12]: P. L. Walters, T. C. Allen, and N. Makri, “Direct determination of discrete harmonic bath parameters from molecular dynamics simulations”, J. Comput. Chem. 38, 110–115 (2017).
[^13]: Z. Tang, X. Ouyang, Z. Gong, H. Wang and J. Wu, "Extended hierarchy equation of motion for the spin-boson model", J. Chem. Phys. 143, 224112 (2015).
[^14]:H. Wang, X. Song, D. Chandler, and W. H. Miller, “Semiclassical study of electronically nonadiabatic dynamics in the condensed-phase: spin-boson problem with Debye spectral density”, J. Chem. Phys. 110, 4828–4840 (1999).
[^15]: M. Thoss, H. Wang, and W. H. Miller, “Self-consistent hybrid approach for complex systems: application to the spin-boson model with Debye spectral density”, J. Chem. Phys. 115, 2991–3005 (2001).
[^16]: I. R. Craig, M. Thoss, and H. Wang, “Proton transfer reactions in model condensed phase environments: accurate quantum dynamics using the multilayer multiconfiguration time-dependent Hartree approach”, J. Chem. Phys. 127, 144503 (2007).
[^17]: A. Ishizaki and G. R. Fleming, "Theoretical examination of quantum coherence in a photosynthetic system at physiological temperature", Proc. Natl. Acad. Sci. U. S. A. 106, 17255-17260 (2009).
[^18]:W.-L. Chan, T. C. Berkelbach, M. R. Provorse, N. R. Monahan, J. R. Tritsch, M. S. Hybertsen, D. R. Reichman, J. Gao and X.-Y. Zhu, "The quantum coherent mechanicsm for singlet fission: Experiment and theory", Acc. Chem. Res. 46, 1321-1329 (2013).
[^19]:G. H. Tao, "Electronically nonadiabatic dynamics in singlet fission: A quasi-classical trajectory simulation", J. Phys. Chem. C 118, 17299-17305 (2014).
[^20]: J. J. Ren, W. T. Li, T. Jiang, Y. H. Wang and Z. G. Shuai, "Time-dependent density matrix renormalization group method for quantum dynamics in complex systems", Wiley Interdiscip. Rev. Comput. Mol. Sci. 12, e1614 (2022).
[^21]: Y. Tanimura and R. Kubo, "Time evolution of a quantum system in contact with a nearly Gaussian-Markoffian noise bath", J. Phys. Soc. Jpn. 58, 101-114 (1989).

[^22]: N. M. Hoffmann, C. Schäfer, A. Rubio, A. Kelly, H. Appel, "Capturing vacuum fluctuations and photon correlations in cavity quantum electrodynamics with multitrajectory Ehrenfest dynamics.", Phys. Rev. A 99, 063819 (2019). 
[^23]: N. M. Hoffmann, C. Schäfer,N. Säkkinen, A. Rubio, H. Appel, A. Kelly, "Benchmarking semiclassical and perturbative methods for real-time simulations of cavity-bound emission and interference", J. Chem. Phys. 151, 244113 (2019).

[^24]: T. Holstein, “Studies of polaron motion: Part II. The “small” polaron”, Ann. Phys. 8, 343–389 (1959).
[^25]: K. Hannewald and P. Bobbert, “Ab initio theory of charge-carrier conduction in ultrapure organic crystals”, Appl. Phys. Lett. 85, 1535–1537 (2004).
[^26]: J. Liu, “Isomorphism between the multi-state Hamiltonian and the second-quantized many-electron Hamiltonian with only 1-electron interactions”, J. Chem. Phys. 146, 024110 (2017).
[^27]: W. Li, J. Ren, and Z. Shuai, “Finite-temperature TD-DMRG for the carrier mobility of organic semiconductors”, J. Phys. Chem. Lett. 11, 4930–4936 (2020).
[^28]: W. Li, J. Ren, and Z. Shuai, “A general charge transport picture for organic semiconductors with nonlocal electron-phonon couplings”, Nat. Commun. 12, 4260 (2021).
[^30]:  X. He, Z. Gong, B. Wu and J. Liu, "Negative zero-point-energy parameter in the Meyer–Miller mapping model for nonadiabatic dynamics", J. Phys. Chem. Lett. 12, 2496–2501 (2021).
[^31]: M. N. Yang and G. R. Fleming, "Influence of phonons on exciton transfer dynamics: Comparison of the Redfield, Forster and modified Redfield equations", Chem. Phys. 275, 355-372 (2002).
[^a]:B. Wu, X. He and J. Liu, "Nonadiabatic Field on Quantum Phase Space: A Century after Ehrenfest", J. Phys. Chem. Lett., 15(2), 644-658 (2024).

[^b]: X. He, X. Cheng, B. Wu and J. Liu, "Nonadiabatic Field with Triangle Window Functions on Quantum Phase Space", J. Phys. Chem. Lett., 15(20), 5452-5466 (2024).

