\page flamelet Laminar flamelets

The flamelet formulation of Peters \cite Peters_1984 is implemented as an option for modeling laminar nonpremixed diffuion flames. This formulation solves the transport equations in the mixture fraction coordinate \f$\xi\f$. The physical configuration is that of an opposed jet with oxidizer in one stream and fuel in the other. Solution in the mixture fraction coordinate assumes species and enthalpy diffuse with equal diffuivities and Lewis numbers of unity,\f$Le_k=\alpha/D_i=k\f$, where \f$\alpha\f$ is the thermal diffusivity and \f$D_k\f$ is the species diffusivity. More generic tranport property relations have been developed \cite Pitsch_1998.

## Species equations

The equations solved for gas species mass fractions \f$y_k\f$ are

$$\prtl{y_k}{t} = \frac{\chi}{2}\prtl{^2y_k}{\xi^2} + \frac{\doot{m}_k^\tprime}{\rho},$$

Where \f$\rho\f$ is density, \f$\chi\f$ is the scalar dissipation rate, and \f$\doot{m}_k^\tprime\f$ is the species reaction rate per unit volume.

The scalar dissipation rate \f$\chi\f$ relates the spatial and mixture fraction coordinates. It is given by
$$\chi(\xi) = 2D\left(\prtl{\xi}{x}\right)^2.$$ The disipation rate profile appearing in the species equation is specified consistent with an opposed jet flame,
$$\chi = \chi_0\exp(-2(\text{erf}^{-1}(2\xi-1))^2).$$
A good approximation to this profile that is used here is given by
$$\chi = \chi_0(1-(2\xi-1)^2)^2,$$
which follows from a hyperbolic tangent profile of \f$\xi(x)\f$, varying from zero at low \f$x\f$ to unity at high \f$x\f$. The prefactor \f$\chi_0\f$ sets the magnitude of the profile and is the value of \f$\chi(\xi=0.5)\f$.

## Energy equation

For adiabatic flames, the enthalpy profile is simply linear between the two streams:
$$h(\xi) = (1-\xi)h_{\xi=0} + (\xi)h_{\xi=1}.$$ This can be generalized for nonadiabatic flames with a radiative source term.

The equation in terms of temperature is given by 
$$\prtl{T}{t} = \frac{\chi}{2}\prtl{^2T}{\xi^2} - \frac{1}{\rho c_p}\sum_kh_k\doot{m}_k^\tprime + \left[\frac{\chi}{2c_p}\prtl{T}{\xi}\prtl{c_p}{\xi} + \frac{\chi}{2c_p}\sum_k\prtl{y_k}{\xi}\prtl{h_k}{\xi}\right] + \frac{Q_r}{\rho c_p}.$$
Here, \f$Q_r\f$ is the radiative source term, discussed in the [diffusion flame](@ref diffusion) section.
This equation was implemented and verified to give the expected linear enthalpy profile when \f$Q_r=0\f$.

## Soot equations

The assumption of a unity Lewis number for soot is not appropriate as soot has negligible diffusion by Brownian motion. It does diffuse by thermophoresis as noted in the [diffusion flame](@ref diffusion) section. A flamelet equation allowing for soot transport by thermophoresis is given by
$$\prtl{m_k}{t} = C\prtl{m_k}{\xi} + S_k.$$
Here, \f$m_k=M_k/\rho\f$, where \f$M_k\f$ is the \f$k^\text{th}\f$ mass-moment. The factor \f$C\f$ can be interpreted as an advection velocity in the mixture fraction coordinate, and is given by
$$C = \frac{0.556\mu}{\rho T}\beta^2\prtl{T}{\xi} - \frac{\beta}{\rho}\prtl{\rho D\beta}{\xi}.$$
\f$\beta\f$ is defined to be
$$\beta\equiv\prtl{\xi}{x} = \sqrt{\frac{\chi}{2D}}.$$
The soot gradient is computed with an upwind discretization using the local value of \f$C\f$.
The source term \f$S_k\f$ is given by
$$S_k = \frac{0.556 \beta m_k}{\rho}\left[\prtl{T}{\xi}\prtl{\mu\beta/T}{\xi} + \frac{\beta\mu}{T}\prtl{^2T}{\xi^2}\right] + \frac{\doot{M}_k}{\rho}.$$

See \cite Lignell_thesis for details and derivations of these equations.

