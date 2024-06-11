#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include "streams.h"
#include "rad_planck_mean.h"
#include "linearInterp.h"
#include "sootHeaders.h"

#include <highfive/highfive.hpp>

#include <vector>
#include <memory>
#include <string>

///////////////////////////////////////////////////////////////////////////////
///
/// One-dimensional flame solver: burner stabilized premixed and diffusion flames
/// Premixed flames can use a fixed temperature profile, or solve the energy equation.
/// Diffusion flames are diffusion only, with Dirichlet boundary conditions.
///   Flame stain occurs through the domain length.
///
///////////////////////////////////////////////////////////////////////////////

class ignis {

public:

    size_t ngrd;                                    ///< number of interior grid points
    size_t nsp;                                     ///< number of gas species
    size_t nsoot;                                   ///< number of soot variables
    size_t nvar;                                    ///< number of transported variables at each grid point
    size_t nvarA;                                   ///< nvar * ngrd

    double                            P;            ///< system pressure, uniform (Pa)
    std::vector<std::vector<double> > y;            ///< mass fractions: y[igrid][isp] 
    std::vector<double>               T;            ///< temperature (K)
    std::vector<std::vector<double> > sootvars;     ///< soot moments or sections; sootvars[igrid][isoot]

    double                            Pstore;       ///< stored system pressure (for initializing from stored state)
    std::vector<std::vector<double> > ystore;       ///< stored mass fractions
    std::vector<double>               Tstore;       ///< stored temperature
    std::vector<std::vector<double> > sootstore;    ///< stored soot variables

    std::vector<double> yLbc, yRbc;                 ///< y boundary values: left and right (as needed)
    double TLbc, TRbc;                              ///< T boundary values: left and right (as needed)
    double hLbc, hRbc;                              ///< h boundary values: left and right (as needed)
    double cpLbc,cpRbc;                             ///< cp boundary values: left and right (as needed)
    std::vector<double> hspLbc;                     ///< species enthalpies on left boundary
    std::vector<double> hspRbc;                     ///< species enthalpies on right boundary

    double              L;                          ///< domain size (m)
    std::vector<double> x;                          ///< grid position values (m)
    std::vector<double> dx;                         ///< grid spacing (m), nonuniform is fine
    std::vector<double> fl;                         ///< fractions for interpolation
    std::vector<double> fr;                         ///< fractions for interpolation

    double Tscale;                                  ///< scaling value for temperature (for solvers)
    double hscale;                                  ///< scaling value for enthalpy (for solvers)
    std::vector<double> sootScales;                 ///< scaling value for soot variables (for solvers)

    std::vector<double> vars0;                      ///< for homotopy approaches
    std::vector<double> F0;                         ///< for homotopy approaches
    double s;                                       ///< homotopy variable

    std::shared_ptr<Cantera::ThermoPhase> gas;      ///< Cantera thermo object
    std::shared_ptr<Cantera::Kinetics>    kin;      ///< Cantera kinetics object
    std::shared_ptr<Cantera::Transport>   trn;      ///< Cantera transport object

    std::shared_ptr<streams> strm;                  ///< streams object (mixture fraction, etc.)
    std::shared_ptr<rad> radProps;                  ///< radiation object
    bool doRadiation;                               ///< radiation flag

    bool doLe1 = false;                             ///< true if doing unity Lewis numbers (default false)  

    double Ttarget;                                 ///< for unsteady cases, run until this max T instead of for a given time
    double dT;                                      ///< delta T increment for unsteady cases
    int isave;                                      ///< file counter for save during unsteady cases

    std::vector<std::vector<double> > flux_y;       ///< species fluxes: [I(igrid, ksp)]    I(igrid,ksp) maps 2D onto 1D
    std::vector<std::vector<double> > flux_soot;    ///< species fluxes: [I(igrid, ksoot)]
    std::vector<double>               flux_h;       ///< species fluxes: [igrid]

    bool isFlamelet = false;                        ///< true for laminar flamelet (mixture fraction coordinate)
    std::vector<double> chi;                        ///< dissipation rate profile
    double chi0;

    bool   isPremixed = false;                      ///< true of the case is a premixed flame, (only left boundary condition, constant mass flux through domain)
    double mflux = 0.0;                             ///< premixed flame mass flux (kg/m2*s)

    bool doEnergyEqn = true;                        ///< for premixed flames: can solve energy equation or set T profile
    std::shared_ptr<linearInterp> LI;               ///< interpolator for specified temperature profiles
    std::vector<double> Tprof_h;                    ///< temperature profile position (h is height above burner (m))
    std::vector<double> Tprof_T;                    ///< temperature profile T values

    //---------------------

    bool doSoot = false;                            ///< soot flag
    std::shared_ptr<soot::sootModel> SM;            ///< soot model
    std::shared_ptr<soot::state>     SMstate;       ///< holds state variables (gas and soot) for soot model

    //---------------------

    std::shared_ptr<HighFive::File> fh5;            ///< hdf5 file pointer

    ////////////////////// member functions

    void setIC(const std::string icType, const std::string fname="");
    void storeState();
    void setFluxesUnity();
    void setFluxes();
    void setGrid(double _L);
    void writeFileHdf5(const std::string gname, const std::string timeType);
    void writeFile(const std::string fname);
    void solveSS();
    void setChi(const double _chi0);
    void solveUnsteady(const double nTauRun, const int nsteps, const bool doWriteTime=true, 
                       const double Tmin=0, const double Tmax=0);
    int  Func(const double *vars, double *F);
    int  rhsf(const double *vars, double *dvarsdt);
    int  rhsf_flamelet(const double *vars, double *dvarsdt);
    void setQrad(std::vector<double> &Q);
    void setTprof(const std::vector<double> &_Tprof_h, const std::vector<double> &_Tprof_T) {
        Tprof_h = _Tprof_h;
        Tprof_T = _Tprof_T;
        LI = std::make_shared<linearInterp>(Tprof_h, Tprof_T);
    }
    void setDerivative2(const double vL, const double vR, 
                        const std::vector<double> &v, 
                        std::vector<double> &d2vdx2);
    void setDerivative( const double vL, const double vR, 
                        const std::vector<double> &v, 
                        std::vector<double> &dvdx);

    size_t I( size_t i, size_t k) { return i*nsp  + k; }     // y[I(i,k)] in 1D --> y[i,k] in 2D
    size_t Ia(size_t i, size_t k) { return i*nvar + k; }     // for indexing combined (a for all) vars

    ////////////////////// constructors 

    ignis(const bool _isPremixed, const bool _doEnergyEqn, const bool _isFlamelet, const bool _doSoot, 
          const size_t _ngrd, const double _L, const double _P,
          std::shared_ptr<Cantera::Solution> csol,
          const std::vector<double> &_yLbc, const std::vector<double> &_yRbc, 
          const double _TLbc, const double _TRbc,
          std::shared_ptr<soot::sootModel> _SM, 
          std::shared_ptr<soot::state>     _SMstate);

};
