#include "flame.h"
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
/// Driver for a simple diffusion flame.
///
///////////////////////////////////////////////////////////////////////////////

int driver_diffusion() {
    
    auto csol = Cantera::newSolution("gri30.yaml");
    auto gas  = csol->thermo();

    //===================== read input file

    YAML::Node inputFile = YAML::LoadFile("../input/input_diffusion.yaml");

    //---------------------

    bool isPremixed = inputFile["isPremixed"].as<bool>();

    //---------------------

    size_t ngrd    = inputFile["ngrd"].as<size_t>();
    double L       = inputFile["L"].as<double>();
    double nTauSS  = inputFile["nTauSS"].as<double>();
    int    nsaveSS = inputFile["nsaveSS"].as<int>();

    vector<double> Ls;
    for(size_t i=0; i<inputFile["Ls"].size(); i++)
        Ls.push_back(inputFile["Ls"][i].as<double>());

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
        SM->coag->set_FM_multiplier(9.0/2.0/2.2);
        //vector<double> sootScales(nsoot, 1.0);
        //sootScales[0] = 1e16;
        //sootScales[1] = 0.01;
        // need a scaling factor for each soot moment
        SMstate = make_shared<state>(nsoot);
        //SMstate->setSootScales(sootScales);
    }

    //---------------------

    bool doEnergyEqn = true;

    //=====================

    flame flm(isPremixed, doEnergyEqn, doSoot, 
              ngrd, L, P, csol,
              yLbc, yRbc, TLbc, TRbc,
              SM, SMstate);

    //flm.doLe1 = true;

    flm.setIC("equilibrium");
    flm.writeFile("IC.dat");

    flm.doRadiation = false;
    flm.solveUnsteady(nTauSS, nsaveSS, false);

    stringstream ss; ss << "L_" << L << "S_" << setfill('0') << setw(3) << 0 << ".dat";
    string fname = ss.str();
    flm.writeFile(fname);

    //---------------

    return 0;
}
