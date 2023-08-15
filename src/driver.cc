#include "flame.h"
#include "cantera/base/Solution.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

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

    flm.writeFile("IC.dat");

    //---------------

    // flm.setGrid(0.03); cout << endl << "L = " << flm.L << endl;
    // flm.solveUnsteady();
    // flm.writeFile("L_0.03U.dat");

    // flm.setGrid(0.02); cout << endl << "L = " << flm.L << endl;
    // flm.solveUnsteady();
    // flm.writeFile("L_0.02U.dat");

    vector<double> Ls = {0.2, 0.04, 0.02, 0.008, 0.006, 0.004, 0.002, 0.001};

    for(int i=0; i<Ls.size(); i++) {
        flm.setGrid(Ls[i]); cout << endl << "L = " << flm.L << endl;
        int nsave = (Ls[i] != 0.002) ? 1 : 40;
        int ntaurun = (Ls[i] != 0.002) ? 5 : 2;
        flm.solveUnsteady(ntaurun, nsave);
        // stringstream ss; ss << "L_" << Ls[i] << "U.dat";
        // string fname = ss.str();
        // flm.writeFile(fname);
    }

    // double Lmax = 0.2;
    // double Lmin = 0.001;
    // double dL = 0.001;
    // for(double L=Lmax; L>=Lmin; L-= dL){
    //     flm.setGrid(L); cout << endl << "L = " << flm.L << endl;
    //     flm.solveSS();
    //     stringstream ss; ss << "L_" << L << ".dat";
    //     string fname = ss.str();
    //     flm.writeFile(fname);
    //}

    // flm.setIC("equilibrium");

    // flm.setGrid(0.2); cout << endl << "L = " << flm.L << endl;
    // flm.solveSS();
    // flm.writeFile("L_0.2.dat");



    // flm.setGrid(0.001); cout << endl << "L = " << flm.L << endl;
    // flm.solveSS();
    // flm.writeFile("L_0.001.dat");



    // double Lmax = 0.2;
    // double Lmin = 0.003;
    // double dL = 0.001;
    // for(double L=Lmax; L>=Lmin; L-= dL){
    //     flm.setGrid(L); cout << endl << "L = " << flm.L << endl;
    //     flm.solveSS();
    //     stringstream ss; ss << "L_" << L << ".dat";
    //     string fname = ss.str();
    //     flm.writeFile(fname);
    // }


    return 0;
}
