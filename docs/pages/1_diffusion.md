\page diffusion Diffusion flames

Diffusion flames are modeled on a simple one-dimensional domain of fixed size \f$L\f$. Only diffusive transport is considered, with no advective transport. The system pressure is assumed spatially uniform and constant in time. By convention, the lower boundary is taken to be oxidizer and the upper boundary is taken to be fuel. This results in the mixture fraction increasing monotonically from zero to one from the lower to the upper boundary. The mixture fraction is the local mass fraction of gas material originating in the fuel stream. 

As noted, there is no advective flame strain, but the fixed domain size imposes a diffusive strain, and at steady state there is a balance between reaction and diffusion. As the domain size is decreased, the diffusive mixing rate increases (through the larger imposed gradients), and the flame eventually extinguishes when the reaction rates cannot keep up with the mixing rates.

This flame configuration was advocated by Pierce for application as a subgrid model for turbulent combustion \cite Pierce_2004, as it is neutral regarding opposed or anti-opposed flow, in contrast to, e.g., [laminar flamelets](@ref flamelet) that assume opposed flow.

## Species equations

The equations solved for gas species mass fractions \f$y_k\f$ are

$$\prtl{y_k}{t} = -\frac{1}{\rho}\prtl{j_k}{x} + \frac{\doot{m}_k^\tprime}{\rho},$$

Where \f$\rho\f$ is density, \f$j_k\f$ is the mass diffusion flux, and \f$\doot{m}_k^\tprime\f$ is the species reaction rate per unit volume.

The unsteady term on the left-hand side of the species equation is retained to allow unsteady evolution to a desired steady state, that is, for computational convenience. As written, the equation neglects the time dependency of density. At steady state, the obtained solution is correct for the given configuration, so the time-independent assumption on density only affects the intermediate unsteady profiles.

<!--
The primary purpose of the unsteady flame is for generation of lookup tables for application to turbulent flows in which the diffusion flames are an idealization of a subrid, for which details of the unsteady evolution not crucial, and the mapping to appropriate state variables is of primary interest. Of course, the assumption can be relaxed, but then the unsteady neglect of advection terms should also be considered.
-->

Dirichlet boundary conditions are specified for species, that is, the species mass fractions are specified.

Ideal gases are assumed, and the density is computed as \f$\rho=MP/RT\f$, where \f$M\f$ is the mean molecular weight of the gas, \f$T\f$ is temperature, \f$P\f$ is pressure, and \f$R\f$ is the gas constant. The diffusion flux \f$j_k\f$ is given by
$$j_k = -\rho D_k\prtl{y_k}{x} - \left(\frac{\rho D_ky_k}{M}\prtl{M}{x}\right),$$
where \f$D_k\f$ is the effective species diffusivity. Mixture average transport properties are assumed. An option for unity Lewis numbers is available, in which case the species diffusivities are all equal to the thermal diffusivity, and the term in parentheses for \f$j_k\f$ is ignored. 

## Energy equation

Energy can be solved in several forms, usually in terms of either enthalpy or temperature. The equation is simpler when solving for enthalpy, but this requires regular inversion to compute the temperature for evaluation of themochemical properties. This inversion requires solution of a nonlinear equation, converged to some tolerance. This is trivial, but when combined with coupled spatial solvers on stiff systems, such as done here, using, e.g., CVODE, the system is less stable than solving temperature directly. Here, the temperature equation is solved, and is given by
$$\prtl{T}{t} = -\frac{1}{\rho c_p}\prtl{q}{x} - \frac{1}{c_p}\sum_kh_k\prtl{y_k}{t} + \frac{Q_r}{\rho c_p},$$
where \f$q\f$ is the heat flux, \f$h_k\f$ is the species enthalpy, \f$Q_r\f$ is the radiative source term (if radiation is used), and \f$c_p\f$ is heat capacity.

This equation for temperature is equivalent to the enthalpy equation given by 
$$\prtl{h}{t} = -\frac{1}{\rho}\prtl{q}{x} + \frac{Q_r}{\rho}.$$
This equivalence was verified in the code implementation by direct comparison.

The heat flux \f$q\f$ is given by
$$q = -k\prtl{T}{x} + \sum_kh_kj_k.$$

[Cantera](https://cantera.org/) is used for all thermochemical and transport properties, including diffusivities, thermal conductivity, viscosity, heat capacity, density, and chemical reaction rates, using available detailed kinetic mechanisms. SI units are used for all quantities (except kmoles are used instead of moles).

### Radiation
The radiative source term \f$Q_r\f$ is computed using an optically-thin approximation,
$$Q_r = -4\sigma k_a(T^4-T_0^4),$$
where \f$\sigma\f$ is the Stefan-Boltzmann constant, \f$k_a\f$ is the absorption coefficient, and \f$T_0\f$ is the temperature of the lower boundary, corresponding nominally to the ambient surrounding air. The temperature and composition dependent absorption coefficient is computed by [RadLib](https://github.com/BYUignite/RadLib.git) \cite Stephens_2022. By default, a Planck-mean assumption is used considering species H\f$_2\f$O, CO\f$_2\f$, CO, and CH\f$_4\f$. Spectral models can also be easily used, including the weighted sum of gray gases (WSSG), and the rank correlated spectral line weighted sum of gray gases (RCSLW).

## Soot

The equation for soot transport is give by

$$\prtl{M_k}{t} = -\prtl{j_{M,k}}{x} + {S_{M_k}},$$

Where \f$M_k\f$ is the \f$k^{th}\f$ soot mass-moment. Alternatively, it can be defined as the soot section \f$k\f$ of size \f$M_k\f$ and corresponding to the number of soot particles per \f$\text{m}^3\f$. \f${S_{M_k}}\f$ is the soot source term and is taken directly from [SootLib](https://github.com/BYUignite/sootlib) \cite Stephens_2023. It includes nucleation, growth, oxidation, and coagulation. \f$M_k\f$ is defined as

$$ M_k = \int m^k n(m) dm $$

where \f$n(m)\f$ is the number of soot particles per \f$\text{m}^3\f$ per kg of soot. Additionally, \f$m\f$ is the mass in kg per soot particle. The value of \f$M_k\f$ when \f$ k = 0 \f$ is taken to be the number of soot particles per \f$\text{m}^3\f$. Likewise, when \f$ k = 1 \f$, \f$ M_k \f$ is representative of the kg of soot per \f$\text{m}^3\f$, which is equal to the denisty of the soot times the mass fraction of the soot. Further values of \f$M_k\f$ can be shown, but do not have as interesting of physical interpretations. The average size of a soot particle, \f$\langle m \rangle \f$ is defined as \f$\frac{M_1}{M_0}\f$, or the average mass in kg per soot particle. Each soot mass-moment is spaced from its predecessor by about 20 orders of magnitude. Thus, to help solve this equation, the soot moments are multiplied by a scaling factor so their values are all about the same order of magnitude. This helps the solver converge, and after the solve, the mass moments are divided by the same scaling factor to revert them back to their original order of magnitude.

The flux of each soot mass moment, \f$j_{M_{k}}\f$ is given by

$$ j_{M_{k}} = -0.556 \nu M_k \frac{\nabla T}{T}, $$

where \f$\nu\f$ is the kinematic viscosity of the gas. This makes the soot transport equation a hyperbolic PDE, and therefore we solve it by upwinding.

## Numerical solution

The partial differential equations given above are solved using a finite volume (FV) formulation on a nonuniform spatial grid. The species equation at grid cell \f$i\f$ is given by 
$$\prtl{y_{k,i}}{t} = -\frac{j_{k,i}^e-j_{k,i}^w}{\rho_i\Delta x_i} + \frac{\doot{m}_{k,i}^\tprime}{\rho_i},$$
where the superscripts \f$e\f$ and \f$w\f$ denote east and west faces of the grid cell, respectively.

The FV energy equation is given by 
$$\prtl{T_i}{t} = -\frac{q_i^e-q_i^w}{\rho_i c_{p,i}\Delta x_i} - \frac{1}{c_{p,i}}\sum_kh_{k,i}\prtl{y_{k,i}}{t} + \frac{Q_{r,i}}{\rho_i c_{p,i}}.$$

## Default configuration

The default configuration corresponds to the composition of the [TNF piloted nonpremixed jet flames](https://tnfworkshop.org/data-archives/pilotedjet/) (flames C, D, E, and F). The case is set up to use the [GRI3.0](http://combustion.berkeley.edu/gri-mech/version30/text30.html) gas mechanism, which consists of 53 gas species and 325 reactions. The domain size is 4 cm, but this is meant to be varied.
