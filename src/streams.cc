#include "streams.h"
#include <cassert>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
///
/// Constructor
///
///////////////////////////////////////////////////////////////////////////////

streams::streams(shared_ptr<Cantera::Solution> csol,
                 const double _P, 
                 const double _h0, 
                 const double _h1,
                 const vector<double> &_y0, 
                 const vector<double> &_y1) :
    P(_P),
    h0(_h0),
    h1(_h1),
    y0(_y0),
    y1(_y1) {

    //----------

    gas  = csol->thermo(); 
    nspc = gas->nSpecies();

    gCHON.resize(4);
    gCHON[0] =  2.0/gas->atomicWeight(gas->elementIndex("C"));
    gCHON[1] =  0.5/gas->atomicWeight(gas->elementIndex("H"));
    gCHON[2] = -1.0/gas->atomicWeight(gas->elementIndex("O"));
    gCHON[3] = 0.0;

    //----------------

    setStoicMixf();

    getMixtureFraction(&y0[0], true);    // set beta0 and beta1
}

///////////////////////////////////////////////////////////////////////////////
/// Computes the temperature, enthalpy, and composition of mixing among streams.
///
/// @param mixf \input mixture fraction, defines elemental composition.
/// @param ymix \output mass fractions of products of complete combustion.
/// @param hmix \output enthalpy of products of complete combustion.
/// @param Tmix \output temperature of products of complete combustion.
///
///////////////////////////////////////////////////////////////////////////////

void streams::getMixingState(const double mixf, vector<double> &ymix,
                             double &hmix, double &Tmix) {

    hmix = h1*mixf + h0*(1.0-mixf);
    for(int k=0; k<nspc; k++)
        ymix[k] = y1[k]*mixf + y0[k]*(1.0-mixf);

    gas->setMassFractions( &ymix[0] );
    gas->setState_HP(hmix, P, 1.E-10);    // get temperature as Tadiabatic

    Tmix = gas->temperature();
}

///////////////////////////////////////////////////////////////////////////////
/// Computes the temperature, enthalpy, and composition of equilibrium at the given mixf.
///
/// @param mixf \input mixture fraction, defines elemental composition.
/// @param yeq  \output mass fractions
/// @param heq  \output enthalpy
/// @param Teq  \output temperature
///
///////////////////////////////////////////////////////////////////////////////

void streams::getEquilibrium_HP(const double mixf, vector<double> &yeq,
                                      double &heq, double &Teq) {

    //---------- Compute the mixing mass fractions and enthalpy

    yeq.resize(nspc);
    for(int k=0; k<nspc; k++)
        yeq[k] = y1[k]*mixf + y0[k]*(1.0-mixf);
    heq = h1*mixf + h0*(1.0-mixf);

    gas->setMassFractions(&yeq[0]);
    gas->setState_HP(heq, P, 1.E-10);

    gas->equilibrate("HP");
    gas->getMassFractions(&yeq[0]);

    Teq = gas->temperature();

}

///////////////////////////////////////////////////////////////////////////////
/// Computes the enthalpy, and composition of equilibrium at the given mixf for given T
///
/// @param mixf \input mixture fraction, defines elemental composition.
/// @param Teq  \input temperature
/// @param yeq  \output mass fractions
/// @param heq  \output enthalpy
///
///////////////////////////////////////////////////////////////////////////////

void streams::getEquilibrium_TP(const double mixf, double Teq, 
                                      vector<double> &yeq, double &heq ) {

    //---------- Compute the mixing mass fractions and enthalpy

    yeq.resize(nspc);
    for(int k=0; k<nspc; k++)
        yeq[k] = y1[k]*mixf + y0[k]*(1.0-mixf);

    gas->setState_TPY(Teq, P, &yeq[0]);

    gas->equilibrate("TP");
    gas->getMassFractions(&yeq[0]);

    heq = gas->enthalpy_mass();

}

///////////////////////////////////////////////////////////////////////////////
/// Computes the temperature, enthalpy, and composition of complete combustion at the given mixf.
/// For nonpremixed flames (don't do anything funny, like have oxygen in the fuel stream)
///
/// @param mixf \input mixture fraction, defines elemental composition.
/// @param ypcc \output mass fractions of products of complete combustion.
/// @param hpcc \output enthalpy of products of complete combustion.
/// @param Tpcc \output temperature of products of complete combustion.
///
///////////////////////////////////////////////////////////////////////////////

void streams::getProdOfCompleteComb(const double mixf, vector<double> &ypcc,
                                    double &hpcc, double &Tpcc) {

    //---------- Compute the mixing mass fractions

    vector<double> yst(nspc, 0.0);
    for(int k=0; k<nspc; k++)
        yst[k] = y0[k]*(1.0-mixfStoic) + y1[k]*mixfStoic;       // holds unburnt mixing of streams

    //--------- Set gas and element indicicies

    gas->setMassFractions(&yst[0]);

    double xC = gas->elementalMoleFraction(gas->elementIndex("C"));
    double xH = gas->elementalMoleFraction(gas->elementIndex("H"));
    double xO = gas->elementalMoleFraction(gas->elementIndex("O"));
    double xN = gas->elementalMoleFraction(gas->elementIndex("N"));

    vector<double> xst(nspc, 0.0);
    xst[gas->speciesIndex("CO2")] = xC;
    xst[gas->speciesIndex("H2O")] = xH/2.0;
    xst[gas->speciesIndex("N2")]  = xN/2.0;

    gas->setMoleFractions(&xst[0]);
    gas->getMassFractions(&yst[0]);                             // holds products of complete combustion

    ypcc.resize(nspc);
    for(size_t k=0; k<nspc; k++) {                              // compute ypcc
        if(mixf < mixfStoic)
            ypcc[k] = y0[k] + (yst[k] - y0[k])*mixf/mixfStoic;
        else
            ypcc[k] = yst[k] + (y1[k] - yst[k])*(mixf-mixfStoic)/(1.0-mixfStoic);
    }

    gas->setMassFractions(&ypcc[0]);
    hpcc = h0*(1.0-mixf) + h1*mixf;
    gas->setState_HP(hpcc, P, 1.E-10);    // get temperature as Tadiabatic
    Tpcc = gas->temperature();
}

///////////////////////////////////////////////////////////////////////////////
///
/// Set the stoichiometric mixture fraction using Bilger's definition */
///
///////////////////////////////////////////////////////////////////////////////

void streams::setStoicMixf() {

    vector<double> x_c(nspc,0.0);

    double mc0, mc1, mo0, mo1, mh0, mh1;

    vector<double> elemMassFrac0 = getElementMassFracs(&y0[0]);
    vector<double> elemMassFrac1 = getElementMassFracs(&y1[0]);

    mc0 = elemMassFrac0[gas->elementIndex("C")]/
        gas->atomicWeight(gas->elementIndex("C"));
    mc1 = elemMassFrac1[gas->elementIndex("C")]/
        gas->atomicWeight(gas->elementIndex("C"));
    mh0 = elemMassFrac0[gas->elementIndex("H")]/
        gas->atomicWeight(gas->elementIndex("H"));
    mh1 = elemMassFrac1[gas->elementIndex("H")]/
        gas->atomicWeight(gas->elementIndex("H"));
    mo0 = elemMassFrac0[gas->elementIndex("O")]/
        gas->atomicWeight(gas->elementIndex("O"));
    mo1 = elemMassFrac1[gas->elementIndex("O")]/
        gas->atomicWeight(gas->elementIndex("O"));

    mixfStoic = (2.0*mc0 + 0.5*mh0 - mo0) /
        (mo1-mo0 + 2.0*(mc0-mc1) + 0.5*(mh0-mh1));

    cout << endl << "# mixfStoic = m_fuel/(m_fuel+m_air) = " << mixfStoic << endl;

}

///////////////////////////////////////////////////////////////////////////////
/// Sets the elements to have the correct Mass Fractions based on the specified array.
///
/// @param y \input mass fraction array to use to get corresponding element fractions.
/// @return vector of element mass fractions.
///
///////////////////////////////////////////////////////////////////////////////

vector<double> streams::getElementMassFracs(const double *y) {

    vector<double> elemMassFrac(gas->nElements(), 0.0);
    gas->setMassFractions( &y[0] );
    for(int k=0; k<gas->nElements(); k++)
        elemMassFrac[k] = gas->elementalMassFraction(k);
    return elemMassFrac;
}

///////////////////////////////////////////////////////////////////////////////
/// Get amount of moles for each element.
///
/// @param x            \input pointer to vector of species mole fractions.
/// @param nOnotFromO2  \input number of moles of oxygen not from O2 (oxygen in the base fuel).
/// @param nHnotFromH2O \input number of moles of hydrogen not from H2O.
/// @param nCnotFromCO2 \input number of moles of carbon not from CO2.
/// @return vector of element moles.
///
///////////////////////////////////////////////////////////////////////////////

vector<double> streams::getElementMoles(const double *x,
                                        double &nOnotFromO2,
                                        double &nHnotFromH2O,
                                        double &nCnotFromCO2) {

    nOnotFromO2  = 0.0;
    nHnotFromH2O = 0.0;
    nCnotFromCO2 = 0.0;

    //vector<double> atomArr(gas->nElements());
    vector<double> elemM(gas->nElements(), 0.0);
    int iO2  = gas->speciesIndex("O2");
    int iO   = gas->elementIndex("O");
    int iCO2 = gas->speciesIndex("CO2");
    int iC   = gas->elementIndex("C");
    int iH2O = gas->speciesIndex("H2O");
    int iH   = gas->elementIndex("H");

    for(int k=0; k<nspc; k++) {
        vector<double> xx(gas->nSpecies(), 0.0);
        xx[k] = 1.0;
        gas->setMoleFractions(&xx[0]);
        for(int m=0; m<gas->nElements(); m++)
            elemM[m] += x[k] * gas->elementalMoleFraction(m);
        if(k != iO2)  nOnotFromO2  += gas->elementalMoleFraction(iO) * x[k];
        if(k != iCO2) nCnotFromCO2 += gas->elementalMoleFraction(iC) * x[k];
        if(k != iH2O) nHnotFromH2O += gas->elementalMoleFraction(iH) * x[k];
    }
    return elemM;
}

///////////////////////////////////////////////////////////////////////////////
/// Compute the mixture fraction from the mass fractions using Bilger's mixf.
/// Set doBeta01=true on first call to initialize members beta0, beta1.
/// Later calls of this function only use the first parameter.
///
/// @param y        \input vector of species mass fractions.
/// @param doBeta01 \input flag=true on first call to set members beta0, beta1.
/// @return mixture fraction
///
///////////////////////////////////////////////////////////////////////////////

double streams::getMixtureFraction(const double *y, const bool doBeta01) {

    vector<double> elemMF;
    double         beta;

    elemMF = getElementMassFracs(y);
    beta   = gCHON[0] * elemMF[gas->elementIndex("C")] +
             gCHON[1] * elemMF[gas->elementIndex("H")] +
             gCHON[2] * elemMF[gas->elementIndex("O")] +
             gCHON[3] * elemMF[gas->elementIndex("N")];

    if(doBeta01) {
        elemMF = getElementMassFracs(&y0[0]);
        beta0  = gCHON[0] * elemMF[gas->elementIndex("C")] +
                 gCHON[1] * elemMF[gas->elementIndex("H")] +
                 gCHON[2] * elemMF[gas->elementIndex("O")] +
                 gCHON[3] * elemMF[gas->elementIndex("N")];

        elemMF = getElementMassFracs(&y1[0]);
        beta1  = gCHON[0] * elemMF[gas->elementIndex("C")] +
                 gCHON[1] * elemMF[gas->elementIndex("H")] +
                 gCHON[2] * elemMF[gas->elementIndex("O")] +
                 gCHON[3] * elemMF[gas->elementIndex("N")];
    }

    if( beta1-beta0 == 0 )
        return 0.0;            // to avoid division by zero
    else
        return (beta - beta0)/(beta1 - beta0);
}
