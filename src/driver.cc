#include "flame.h"
#include "cantera/base/Solution.h"
#include "yaml-cpp/yaml.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int main() {
    
    auto csol = Cantera::newSolution("gri30.yaml");
    auto gas  = csol->thermo();

    //---------------------

    YAML::Node inputFile = YAML::LoadFile("../input.yaml");

    size_t ngrd = inputFile["ngrd"].as<size_t>();
    double L = inputFile["L"].as<double>();
    double nTauRun = inputFile["nTauRun"].as<double>();
    size_t nSteps = inputFile["nSteps"].as<size_t>();

    bool   isPremixed  = inputFile["isPremixed"].as<bool>();
    bool   doRadiation = inputFile["doRadiation"].as<bool>();

    size_t P = inputFile["P"].as<double>();
    double v = inputFile["LBC"]["v"].as<double>();
    double TLbc = inputFile["LBC"]["TLbc"].as<double>();
    vector<double> xLbc(gas->nSpecies());
    YAML::Node xx = inputFile["LBC"]["comp"];
    for(auto it=xx.begin(); it!=xx.end(); it++)
        xLbc[gas->speciesIndex(it->first.as<string>())] = it->second.as<double>();

    bool doEnergyEqn = inputFile["doEnergyEqn"].as<bool>();
    vector<double> Tprof_h;
    vector<double> Tprof_T;
    if(!doEnergyEqn) {
        for(size_t i=0; i<inputFile["Tprof"].size(); i++) {
            Tprof_h.push_back(inputFile["Tprof"][i][0].as<double>());
            Tprof_T.push_back(inputFile["Tprof"][i][1].as<double>());
        }
    }

    //---------------------

    gas->setState_TPX(TLbc, P, &xLbc[0]);
    vector<double> yLbc(ngrd);
    gas->getMassFractions(&yLbc[0]);
    double mflux = gas->density()*v;

    flame flm(isPremixed, doEnergyEqn, ngrd, L, P, csol,
              yLbc, yLbc, TLbc, TLbc);
    flm.mflux = mflux;
    if(!doEnergyEqn) flm.setTprof(Tprof_h, Tprof_T);

    flm.setIC("premixed");
    flm.writeFile("IC.dat");

    //---------------------

    flm.solveUnsteady(nTauRun, nSteps, false);
    string fname = "premixed.dat";
    flm.writeFile(fname);

    //---------------------

    return 0;
}
