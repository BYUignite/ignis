#include "ignis.h"
#include "cantera/base/Solution.h"
#include "yaml-cpp/yaml.h"
#include "sootHeaders.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace soot;

///////////////////////////////////////////////////////////////////////////////
///
/// Driver for a simple laminar flamelet solved in mixture fraction space
///
///////////////////////////////////////////////////////////////////////////////

int driver_flamelet() {
    
    // auto csol = Cantera::newSolution("gri30.yaml");
    //auto csol = Cantera::newSolution("../input/c2h4det.yaml");
    auto csol = Cantera::newSolution("../input/LUsk17.yaml");
    auto gas  = csol->thermo();

    //===================== read input file

    YAML::Node inputFile = YAML::LoadFile("../input/input_flamelet.yaml");

    //---------------------

    bool isFlamelet  = true;
    bool isPremixed  = false;
    bool doEnergyEqn = inputFile["doEnergyEqn"].as<bool>();
    bool doUnifChi   = inputFile["doUnifChi"].as<bool>();

    double L       = 1.0;       // mixture fraction ranges from 0 to 1
    size_t ngrd    = inputFile["ngrd"].as<size_t>();
    double nTauSS  = inputFile["nTauSS"].as<double>();
    int    nsaveSS = inputFile["nsaveSS"].as<int>();

    double chi0    = inputFile["chi0"].as<double>();

    //--------------------- gas streams

    double P = inputFile["P"].as<double>();

    double TLbc = inputFile["LBC"]["TLbc"].as<double>();
    vector<double> yLbc(gas->nSpecies());
    YAML::Node yy = inputFile["LBC"]["comp"];
    for(auto it=yy.begin(); it!=yy.end(); it++)
        yLbc[gas->speciesIndex(it->first.as<string>())] = it->second.as<double>();

    double TRbc = inputFile["RBC"]["TRbc"].as<double>();
    vector<double> yRbc(gas->nSpecies());
    yy = inputFile["RBC"]["comp"];
    for(auto it=yy.begin(); it!=yy.end(); it++)
        yRbc[gas->speciesIndex(it->first.as<string>())] = it->second.as<double>();

    //--------------------- soot

    bool   doSoot = inputFile["doSoot"].as<bool>();
    size_t nsoot  = doSoot ? inputFile["nsoot"].as<size_t>() : 0;

    shared_ptr<sootModel> SM;
    shared_ptr<state>     SMstate;

    if(doSoot) {

        nucleationModel  *nucl = new soot::nucleationModel_LL();
        growthModel      *grow = new soot::growthModel_LL();
        oxidationModel   *oxid = new soot::oxidationModel_LL();
        coagulationModel *coag = new soot::coagulationModel_FM();

        SM = make_shared<sootModel_QMOM>(nsoot, nucl, grow, oxid, coag);
        //SM->coag->set_FM_multiplier(9.0/2.0/2.2);
        SM->coag->set_FM_multiplier(9.0/2.2);
        SMstate = make_shared<state>(nsoot);
    }

    //--------------------- radiation

    string  radType = inputFile["radType"]  ?  inputFile["radType"].as<string>() : "planckmean";
    
    //---------------------


    //=====================

    ignis flm(isPremixed, doEnergyEqn, isFlamelet, doSoot, 
              ngrd, L, P, csol, radType,
              yLbc, yRbc, TLbc, TRbc,
              SM, SMstate);

    flm.doUnifChi = doUnifChi;
    flm.setChi(chi0);

    flm.setIC("equilibrium");
    flm.writeFile("IC.dat");

    flm.doRadiation = false;
    flm.solveUnsteady(nTauSS, nsaveSS, true);

    stringstream ss; ss << "X_" << chi0 << "S_" << setfill('0') << setw(3) << 0 << ".dat";
    string fname = ss.str();
    flm.writeFile(fname);

    //---------------

    return 0;
}
