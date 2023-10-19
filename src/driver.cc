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

    size_t ngrd = 40;
    double L = 0.1;
    double P = 101325;

    double TLbc = 300.0;
    vector<double> yLbc(gas->nSpecies());
    yLbc[gas->speciesIndex("O2")] = 0.233;
    yLbc[gas->speciesIndex("N2")] = 0.767;

    double TRbc = 300.0;
    vector<double> yRbc(gas->nSpecies());
    yRbc[gas->speciesIndex("CH4")] = 1;

    flame flm(ngrd, L, P, csol,
              yLbc, yRbc, TLbc, TRbc);

    flm.setIC("equilibrium");
    flm.storeState();
    flm.writeFile("IC.dat");

    //---------------

    double nTauSS     = 5;
    int    nsaveSS    = 1;
    double nTauU      = 5;
    int    nsaveU     = 10;

    vector<double> Ls = {0.2, 0.04, 0.02, 0.008, 0.006, 0.004, 0.002};
    
    double Tmin, Tmax;

    for(int i=0; i<Ls.size(); i++) {
        L = Ls[i];

        //----- do SS solution

        Tmax = *max_element(flm.T.begin(), flm.T.end());
        flm.setGrid(L); cout << "\n\nL = " << flm.L << endl;
        cout << endl << "do SS"; cout.flush();
        flm.LdoRadiation = false;
        flm.solveUnsteady(nTauSS, nsaveSS, false);

        //----- if blows out, rerun, store as it blows out, else run unsteady with heat loss

        if(*max_element(flm.T.begin(), flm.T.end()) < 1.05*min(TLbc, TRbc)) {
            cout << endl << "Extinction for L=" << flm.L << endl;
            Tmin = *max_element(flm.T.begin(), flm.T.end());
            flm.setIC("stored");
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
            flm.LdoRadiation = true;
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
