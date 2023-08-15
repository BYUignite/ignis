#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include "streams.h"

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

    streams strm;

    std::vector<std::vector<double> > flux_y;      // flux_y[I(igrid, ksp)]
    std::vector<double>               flux_h;      // flux_h[igrid]

    ////////////////////// member functions

    void setIC(std::string icType, std::string fname="");
    void setFluxes();
    void setGrid(double _L);
    void writeFile(std::string fname);
    void solveSS();
    void solveUnsteady(int ntaurun, int nsave);
    int  Func(const double *vars, double *F);
    int  rhsf(const double *vars, double *dvarsdt);

    size_t I( size_t i, size_t k) { return i*nsp  + k; }     // y[I(i,k)] in 1D --> y[i,k] in 2D
    size_t Ia(size_t i, size_t k) { return i*nvar + k; }     // for indexing combined (a for all) vars

    ////////////////////// constructors 

    flame(const size_t _ngrd, const double _L, const double _P,
          std::shared_ptr<Cantera::Solution> csol,
          const std::vector<double> &_yLbc, const std::vector<double> &_yRbc, 
          const double _TLbc, const double _TRbc);
};
