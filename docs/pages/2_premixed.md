\page premixed Premixed flames

The premixed flame configuration solves burner-stabilized premixed flames. This is a common and convenient research configuration. The burner is typically a flat, cooled, usually porous surface through which a uniform fuel/oxidizer mixture is passed. The flow rate is smaller than the laminar flame speed, so the flame propagates towards the burner face, but cannot pass through the cooled porous openings. This stabilizes the flame and provides a well-defined boundary condition for solution. The temperature is fixed at the specified cooled burner value. The total mass flow rate (or mass flux \f$\doot{m}^\dprime\f$) is specified by the upstream mass flow controllers. Since species can diffuse upstream of the burner, the total species flux at the burner face is specified (advective and diffusive), rather than the species mass fractions. At the exit of the domain of size \f$L\f$, Neumann boundary conditions with zero gradients are specified.

## Diffusion fluxes

The equations solved for temperature and gas species are the same as for the [diffusion flame](@ref diffusion), but the species and heat fluxes are computed differently. For the premixed flame, the species flux is given by 
$$j_k = \doot{m}^\dprime y_k-\rho D_k\prtl{y_k}{x} - \left(\frac{\rho D_ky_k}{M}\prtl{M}{x}\right).$$
This is the same as the flux for the diffusion flames, except for the addition of the first term \f$\doot{m}^\dprime y_k\f$, which is the advective flux of species \f$k\f$. At the \f$x=0\f$ boundary, the flux is simply taken as \f$j_k=\doot{m}^\dprime y_{k,0}\f$, where \f$y_{k,0}\f$ is the upstream feed-gas mass fraction of the species.

The heat flux is the same as for the diffusion flames, and is given by
$$q = -k\prtl{T}{x} + \sum_kh_kj_k,$$
where \f$j_k\f$ is given above. The temperature gradient at the \f$x=0\f$ face (burner surface), is computed numerically between the center of the first grid cell and the specified burner temperature \f$T_0\f$.

## Energy equation options

For the premixed flame, the energy equation can be solved as described, including a radiative source term if desired, or the temperature profile can be specified directly and no energy equation solved. In the latter case, a discrete temperature profile is specified in the input file, and then linear interpolation is used to evaluate the temperature on the chosen computational grid. Such temperature profiles are commonly applied from measured experimental data in order to remove difficulty in accurate modeling of the profile. This is especially true when radiative effects and/or soot are present, or when soot modeling is the objective of the simulation, since errors in modeling the temperature have a strong impact on modeling the soot profile.

## Soot equations
Soot is treated as for the [diffusion flames](@ref diffusion).

## Verification

The premixed flame code solutions were verified by comparison to output from Cantera for the same flame configuration, with identical results.

## Default configuration

The default configuration corresponds to the [ISF-4 premixed flame](https://www.adelaide.edu.au/cet/isfworkshop/data-sets/laminar-flames#isf-4-premixed-flames-2-mckenna-burner-stabilised-flames) 2: McKenna burner-stabilized flames, flame a) ethylene/air at an equivalence ratio of \f$\phi=2.34\f$, with the experimental temperature profile used (no energy eqauation).
