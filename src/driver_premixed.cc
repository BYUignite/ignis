#include "flame.h"
#include "cantera/base/Solution.h"
#include "yaml-cpp/yaml.h"
#include "sootHeaders.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>

#include <memory>

using namespace std;
using namespace soot;

///////////////////////////////////////////////////////////////////////////////

int driver_premixed() {
    
    auto csol = Cantera::newSolution("gri30.yaml");
    auto gas  = csol->thermo();

    //===================== read input file

    YAML::Node inputFile = YAML::LoadFile("../input/input_premixed.yaml");

    //---------------------

    bool   isPremixed  = inputFile["isPremixed"].as<bool>();

    //---------------------

    size_t ngrd = inputFile["ngrd"].as<size_t>();
    double L = inputFile["L"].as<double>();
    double nTauRun = inputFile["nTauRun"].as<double>();
    size_t nSteps = inputFile["nSteps"].as<size_t>();

    //---------------------

    bool   doSoot      = inputFile["doSoot"].as<bool>();
    size_t nsoot       = doSoot ? inputFile["nsoot"].as<size_t>() : 0;

    shared_ptr<sootModel> SM;
    shared_ptr<state>     SMstate;

    if(doSoot) {

        nucleationModel  *nucl = new soot::nucleationModel_LL();
        growthModel      *grow = new soot::growthModel_LL();
        oxidationModel   *oxid = new soot::oxidationModel_LL();
        coagulationModel *coag = new soot::coagulationModel_FM();

        SM = make_shared<sootModel_QMOM>(nsoot, nucl, grow, oxid, coag);
        SM->coag->set_FM_multiplier(9.0/2.0/2.2);
        //vector<double> sootScales(nsoot, 1.0);
        //sootScales[0] = 1e16;
        //sootScales[1] = 0.01;
        // need a scaling factor for each soot moment
        SMstate = make_shared<state>(nsoot);
        //SMstate->setSootScales(sootScales);
    }

    //---------------------

    double P = inputFile["P"].as<double>();
    double v = inputFile["v"].as<double>();
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

    //=====================

    gas->setState_TPX(TLbc, P, &xLbc[0]);
    vector<double> yLbc(ngrd);
    gas->getMassFractions(&yLbc[0]);
    double mflux = gas->density()*v;

    flame flm(isPremixed, doEnergyEqn, doSoot, 
              ngrd, L, P, csol, 
              yLbc, yLbc, TLbc, TLbc, 
              SM, SMstate);

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
