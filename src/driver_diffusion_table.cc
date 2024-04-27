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
/// Driver to run diffusion flames for various domain sizes L and times t.
/// Adiabatic steady state files at a given L, and for each of those, run unsteady with radiative loss.
///   (except for cases that blows out due to strain.)
/// These files can be used to make a lookup table.
/// L and t can be mapped to progress variable and enthalpy.
///
///////////////////////////////////////////////////////////////////////////////

int driver_diffusion_table() {
    
    auto csol = Cantera::newSolution("gri30.yaml");
    auto gas  = csol->thermo();

    //===================== read input file

    YAML::Node inputFile = YAML::LoadFile("../input/input_diffusion_table.yaml");

    //---------------------

    bool isPremixed = inputFile["isPremixed"].as<bool>();

    //---------------------

    size_t ngrd    = inputFile["ngrd"].as<size_t>();
    double L       = inputFile["L"].as<double>();
    double nTauSS  = inputFile["nTauSS"].as<double>();
    double nTauU   = inputFile["nTauU"].as<double>();
    int    nsaveSS = inputFile["nsaveSS"].as<int>();
    int    nsaveU  = inputFile["nsaveU"].as<int>();

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
    flm.storeState();
    flm.writeFile("IC.dat");

    //---------------

    double Tmin, Tmax;

    for(int i=0; i<Ls.size(); i++) {
        L = Ls[i];

        //----- do SS solution

        Tmax = *max_element(flm.T.begin(), flm.T.end());
        flm.setGrid(L); cout << "\n\nL = " << flm.L << endl;
        cout << endl << "do SS"; cout.flush();
        flm.doRadiation = false;
        flm.solveUnsteady(nTauSS, nsaveSS, false);

        //----- if blows out, rerun, store as it blows out, else run unsteady with heat loss

        if(*max_element(flm.T.begin(), flm.T.end()) < 1.5*min(TLbc, TRbc)) {
            cout << endl << "Extinction for L=" << flm.L << endl;
            Tmin = *max_element(flm.T.begin(), flm.T.end());
            flm.setIC("stored");
            stringstream ss; ss << "L_" << L << "S_" << setfill('0') << setw(3) << 0 << ".dat";
            string fname = ss.str();
            flm.writeFile(fname);
            flm.solveUnsteady(nTauU,  nsaveU, false, Tmin, Tmax);
        }

        else {
            stringstream ss; ss << "L_" << L << "S_" << setfill('0') << setw(3) << 0 << ".dat";
            string fname = ss.str();
            flm.writeFile(fname);

            Tmax = *max_element(flm.T.begin(), flm.T.end());
            flm.storeState();

            //----- do unsteady with radiation to find Tmin

            cout << endl << "do unsteady to get Tmin";
            flm.doRadiation = true;
            flm.solveUnsteady(nTauU,  1, false);
            Tmin = *max_element(flm.T.begin(), flm.T.end());
            flm.setIC("stored");

            //----- do unsteady with radiation nominally evenly spaced between Tmax and Tmin

            cout << endl << "do unsteady";
            flm.solveUnsteady(nTauU,  nsaveU, false, Tmin, Tmax);
            flm.setIC("stored");
        }
    }

    return 0;
}
