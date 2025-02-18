#include "ignis.h"
#include "cantera/base/Solution.h"
#include "yaml-cpp/yaml.h"
#include "sootHeaders.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace soot;

///////////////////////////////////////////////////////////////////////////////
///
/// Driver to run flamelets for various domain sizes L and times t.
///
///////////////////////////////////////////////////////////////////////////////

int driver_flamelet_table() {
    
    auto csol = Cantera::newSolution("../input/LUsk17.yaml");
    auto gas  = csol->thermo();

    //===================== read input file

    YAML::Node inputFile = YAML::LoadFile("../input/input_flamelet_table.yaml");

    //---------------------

    bool isFlamelet  = true;
    bool isPremixed  = false;
    bool doEnergyEqn = false;                             // h(x) profile is specified directly
    bool doUnifChi   = inputFile["doUnifChi"].as<bool>();

    double L       = 1.0;
    size_t ngrd    = inputFile["ngrd"].as<size_t>();
    double nTauSS  = inputFile["nTauSS"].as<double>();
    double nTauU   = inputFile["nTauU"].as<double>();
    int    nsaveSS = inputFile["nsaveSS"].as<int>();
    int    nsaveU  = inputFile["nsaveU"].as<int>();

    vector<double> chiList;
    for(size_t i=0; i<inputFile["chiList"].size(); i++)
        chiList.push_back(inputFile["chiList"][i].as<double>());

    vector<double> hlList;
    for(size_t i=0; i<inputFile["hlList"].size(); i++)
        hlList.push_back(inputFile["hlList"][i].as<double>());

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

    //---------------------

    string radType = "none";
    bool   doSoot  = inputFile["doSoot"] ? inputFile["doSoot"].as<bool>() : false;
    shared_ptr<sootModel> SM;
    shared_ptr<state>     SMstate;

    //---------------------

    ignis flm(isPremixed, doEnergyEqn, isFlamelet, doSoot, 
              ngrd, L, P, csol, radType,
              yLbc, yRbc, TLbc, TRbc,
              SM, SMstate);

    flm.doUnifChi = doUnifChi;
    flm.doRadiation = false;

    double pvMin, pvMax;
    vector<int> ichiExt(hlList.size(), chiList.size());        // for each hl, record ichi that extinguished

    flm.setIC("equilibrium");
    flm.writeFile("IC.dat");
    flm.storeState(1);

    size_t nsaveU_;

    stringstream ss;
    string fname;

    ///////////////////////

    for(int ichi=0; ichi<chiList.size(); ichi++) {           // loop chi values

        cout << "\n\nChi = " << chiList[ichi];
        flm.setChi(chiList[ichi]);
        pvMax = *max_element(flm.pv.begin(), flm.pv.end());

        //----- do adiabatic solution

        cout << endl << "\t" << "adiabatic"; cout.flush();
        flm.set_h(0.0);
        flm.solveUnsteady(nTauSS, nsaveSS, false);

        //----- if blow out, rerun, store as it blows out

        //if(*max_element(flm.T.begin(), flm.T.end()) < 2.0*min(TLbc, TRbc)) {
        if(flm.pvMaxForFlmltExtHl != -1) {
            cout << "... extinction";
            pvMin = flm.pvMaxForFlmltExtHl; // *max_element(flm.pv.begin(), flm.pv.end());
            cout << endl << "\trun unsteady through extinction:";
            cout << endl << "\tpvMin, pvMax: " << pvMin << " " << pvMax;
            flm.setIC("stored_1");
            //flm.solveUnsteady(nTauU,  nsaveU, false, pvMin, pvMax);
            nsaveU_ = nsaveU + (chiList.size()-ichi);
            flm.solveUnsteady(nTauU,  nsaveU_, false, pvMin, pvMax);
            ichiExt[0] = ichi;
            flm.pvMaxForFlmltExtHl = -1;      // reset value
            break;
        }
        else {
            ss.str(""); ss << "X_" << chiList[ichi] << "_hl_0.dat";
            fname = ss.str();
            flm.writeFile(fname);
            if(ichi==0) {
                flm.set_hsens();
                ofstream ofile("hsens.dat");
                ofile << "# mixf, hsens (J/kg)";
                ofile << endl << 0.0 << "  " << flm.hLbc;
                for(int ii=0; ii<ngrd; ii++)
                    ofile << endl << flm.x[ii] << "  " << flm.hsens[ii];
                ofile << endl << 1.0 << "  " << flm.hRbc;
                ofile.close();
            }
            flm.storeState(1);

            //=========================

            for(int ihl=1; ihl<hlList.size(); ihl++) {     // loop heat loss values for given chi
                cout << endl << "\thl = " << hlList[ihl];
                if(ichi >= ichiExt[ihl]){ 
                    cout << "... skipping (previous extinction)";
                    continue;
                }
                flm.set_h(hlList[ihl]);
                flm.storeState(2);
                pvMax = *max_element(flm.pv.begin(), flm.pv.end());
                flm.solveUnsteady(nTauSS, nsaveSS, false);

                //if(*max_element(flm.T.begin(), flm.T.end()) < 2.0*min(TLbc, TRbc)) { // blowout
                if(flm.pvMaxForFlmltExtHl != -1) {
                    cout << "\t" << "... extinction";
                    pvMin = flm.pvMaxForFlmltExtHl; // *max_element(flm.pv.begin(), flm.pv.end());
                    cout << endl << "\trun unsteady through extinction:";
                    cout << endl << "\tpvMin, pvMax: " << pvMin << " " << pvMax;
                    flm.setIC("stored_2");
                    nsaveU_ = nsaveU + (chiList.size()-ichi);
                    flm.solveUnsteady(nTauU,  nsaveU_, false, pvMin, pvMax);
                    ichiExt[ihl] = ichi;
                    flm.pvMaxForFlmltExtHl = -1;      // reset value

                    flm.setIC("stored_2");
                    //break;
                }
                else {
                    ss.str(""); ss << "X_" << chiList[ichi] << "_hl_" << hlList[ihl] << ".dat";
                    fname = ss.str();
                    flm.writeFile(fname);
                
                }
            }
            flm.setIC("stored_1");
        }

    }

    /////////////////////// do Products of Complete Combustion: all heat losses
    // NOTE: PCC does not give a good progress variable since rich products are CO2 and fuel

    //cout << "\n\nWriting files for equilibrium" << endl;

    //flm.setIC("equilibrium");
    //for(int ihl=0; ihl<hlList.size(); ihl++) {     // loop heat loss values for given chi
    //    flm.set_h(hlList[ihl]);
    //    flm.set_T();
    //    ss.str(""); ss << "EQ_hl_" << hlList[ihl] << ".dat";
    //    fname = ss.str();
    //    flm.writeFile(fname);
    //}

    ///////////////////////

    return 0;
}
