#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include "streams.h"
#include "rad_planck_mean.h"

#include <vector>
#include <memory>
#include <string>

///////////////////////////////////////////////////////////////////////////////

class flame {

public:

    size_t ngrd;
    size_t nsp;
    size_t nvar;
    size_t nvarA;

    double                            P;
    std::vector<std::vector<double> > y;           // y[igrid][isp]
    std::vector<double>               T;

    double                            Pstore;
    std::vector<std::vector<double> > ystore;      // y[igrid][isp]
    std::vector<double>               Tstore;

    std::vector<double> yLbc, yRbc;
    double TLbc, TRbc;
    double hLbc, hRbc;

    double              L;
    std::vector<double> x;
    std::vector<double> dx;

    double Tscale;
    double hscale;
    std::vector<double> vars0;
    std::vector<double> F0;
    double s;

    std::shared_ptr<Cantera::ThermoPhase> gas;
    std::shared_ptr<Cantera::Kinetics>    kin;
    std::shared_ptr<Cantera::Transport>   trn;

    std::shared_ptr<streams> strm;
    rad     *planckmean;
    bool    LdoRadiation;

    double Ttarget;
    double dT;
    int isave;

    std::vector<std::vector<double> > flux_y;      // flux_y[I(igrid, ksp)]
    std::vector<double>               flux_h;      // flux_h[igrid]

    bool LisPremixed = false;
    double mflux = 0.0;

    ////////////////////// member functions

    void setIC(std::string icType, std::string fname="");
    void storeState();
    void setFluxesUnity();
    void setFluxes();
    void setGrid(double _L);
    void writeFile(std::string fname);
    void solveSS();
    void solveUnsteady(double nTauRun, int nsteps, bool LwriteTime=true, double Tmin=0, double Tmax=0);
    int  Func(const double *vars, double *F);
    int  rhsf(const double *vars, double *dvarsdt);
    void setQrad(std::vector<double> &Q);

    size_t I( size_t i, size_t k) { return i*nsp  + k; }     // y[I(i,k)] in 1D --> y[i,k] in 2D
    size_t Ia(size_t i, size_t k) { return i*nvar + k; }     // for indexing combined (a for all) vars

    ////////////////////// constructors 

    flame(const bool _LisPremixed, const size_t _ngrd, const double _L, const double _P,
          std::shared_ptr<Cantera::Solution> csol,
          const std::vector<double> &_yLbc, const std::vector<double> &_yRbc, 
          const double _TLbc, const double _TRbc);

    ~flame() {
        delete planckmean;
    }
};
