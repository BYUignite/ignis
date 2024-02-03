#include "flame.h"
#include "cantera/base/Solution.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int main() {
    
    //auto csol = Cantera::newSolution("c2h4red.yaml");
    auto csol = Cantera::newSolution("gri30.yaml");
    auto gas  = csol->thermo();

    size_t ngrd = 100;
    double L = 0.05;
    double P = 101325;

    double TLbc = 298.0;
    vector<double> yLbc(gas->nSpecies());
    yLbc[gas->speciesIndex("C2H4")] = 0.1408;    // phi = 2.34
    yLbc[gas->speciesIndex("O2")]  = 0.1805;
    yLbc[gas->speciesIndex("N2")]  = 0.6787;

    gas->setState_TPY(TLbc, P, &yLbc[0]);
    double rho = gas->density();                  // kg/m3
    double v   = 0.0673;                         // m/s
    double mflux = rho*v;                        // kg/m2*s

    double TRbc = TLbc;
    vector<double> yRbc = yLbc;

    bool LisPremixed = true;

    flame flm(LisPremixed, ngrd, L, P, csol,
              yLbc, yRbc, TLbc, TRbc);

    flm.mflux = mflux;

    flm.setIC("premixed");
    flm.writeFile("IC.dat");

    //---------------

    double nTauRun = 5.0;
    int    nsteps  = 50;

    flm.solveUnsteady(nTauRun, nsteps, false);
    string fname = "premixed.dat";
    flm.writeFile(fname);

    //---------------


    return 0;
}
