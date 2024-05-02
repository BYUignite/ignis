
#include "fuego.h"
#include "cantera/base/ct_defs.h"
#include "solver_kinsol.h"
#include "integrator_cvode.h"

#include <iostream>
#include <algorithm>        // max
#include <iomanip>
#include <fstream>
#include <sstream>
#include <numeric>


using namespace std;
using soot::sootModel, soot::state, soot::gasSp, soot::gasSpMapIS, soot::gasSpMapES;

////////////////////////////////////////////////////////////////////////////////

int Func_kinsol(N_Vector varsKS, N_Vector fvec, void *user_data);
int rhsf_cvode(realtype t, N_Vector varsCV, N_Vector dvarsdtCV, void *user_data);

////////////////////////////////////////////////////////////////////////////////
///
/// Constructor function
/// @param _isPremixed  \input flag indicating if running a premixed or diffusion flame
/// @param _doEnergyEqn \input flag indicating if solving the energy equation or using a given T profile (for premixed flames)
/// @param _doSoot      \input flag indicating if solving with soot
/// @param _ngrd        \input number of interior grid points
/// @param _L           \input domain size (m)
/// @param _P           \input system pressure (Pa)
/// @param csol         \input Cantera Solution object (used to set gas, kin, trn)
/// @param _yLbc        \input mass fractions on left boundary
/// @param _yRbc        \input mass fractions on right boundary (not used for premixed flames)
/// @param _TLbc        \input temperature (K) on left boundary
/// @param _TRbc        \input temperature (K) on right boundary (not used for premixed flames)
/// @param _SM          \input soot model object
/// @param _SMstate     \input soot model state object (gas and soot properties for input to SM)
///
////////////////////////////////////////////////////////////////////////////////

fuego::fuego(const bool _isPremixed, 
             const bool _doEnergyEqn,
             const bool _isFlamelet,
             const bool _doSoot,
             const size_t _ngrd, const double _L, double _P, 
             shared_ptr<Cantera::Solution> csol,
             const vector<double> &_yLbc, const vector<double> &_yRbc, 
             const double _TLbc, const double _TRbc,
             shared_ptr<sootModel> _SM, 
             shared_ptr<state>     _SMstate) :
    isPremixed(_isPremixed),
    isFlamelet(_isFlamelet),
    doEnergyEqn(_doEnergyEqn),
    doSoot(_doSoot),
    ngrd(_ngrd),
    L(_L),
    P(_P),
    yLbc(_yLbc),
    yRbc(_yRbc),
    TLbc(_TLbc),
    TRbc(_TRbc),
    SM(_SM),
    SMstate(_SMstate) {

    //----------

    if(_isPremixed && _isFlamelet)
        throw runtime_error("Cannot set both premixed and flamelet");
    if(_isFlamelet && _L!=1.0)
        throw runtime_error("flamelet should have L=1");

    //----------

    gas = csol->thermo(); 
    kin = csol->kinetics(); 
    trn = csol->transport(); 

    nsp   = gas->nSpecies();
    nsoot = (doSoot ? SM->nsoot : 0);
    nvar  = nsp + nsoot + 1;
    nvarA = nvar*ngrd;

    y = vector<vector<double> >(ngrd, vector<double>(nsp, 0.0));
    if(doSoot) sootvars = vector<vector<double> >(ngrd, vector<double>(nsoot, 0.0));
    T.resize(ngrd, 0.0);

    //---------- set grid

    setGrid(L);

    //----------

    gas->setState_TPY(TLbc, P, &yLbc[0]);
    hLbc = gas->enthalpy_mass();
    cpLbc = gas->cp_mass();
    hspLbc.resize(nsp);
    gas->getEnthalpy_RT(&hspLbc[0]);
    for(size_t k=0; k<nsp; k++)                   // --> hsp = J/kg species i
        hspLbc[k] *= TLbc*Cantera::GasConstant/gas->molecularWeight(k);

    gas->setState_TPY(TRbc, P, &yRbc[0]);
    hRbc = gas->enthalpy_mass();
    cpRbc = gas->cp_mass();
    hspRbc.resize(nsp);
    gas->getEnthalpy_RT(&hspRbc[0]);
    for(size_t k=0; k<nsp; k++)                   // --> hsp = J/kg species i
        hspRbc[k] *= TRbc*Cantera::GasConstant/gas->molecularWeight(k);

    if(!isPremixed) strm = make_shared<streams>(csol, P, hLbc, hRbc, yLbc, yRbc);

    hscale = max(abs(hLbc), abs(hRbc));
    Tscale = 2500;

    //---------- set sootScales, and species map
    // Mk = integral (m^k * n(m) *dm) --> Mk ~ M0<m>^k, where <m> in an average soot size 
    // <m> = M1/M0 = rho*Ys/M0 = fv*rhos/M0
    // Mk ~ M0(rhos*fv/M0)^k; take M0=1E10 #/m3, rhos=1850 kg/m3, fv=1E-6

    if(doSoot) {
        sootScales = vector<double>(nsoot, 1E20);
        for(int i=1; i<nsoot; i++)
            sootScales[i] = sootScales[i-1]*(1850*1E-6/sootScales[0]);
        sootScales = {1E20, 1.8E-3, 3.4E-26, 6.3E-49, 1.2E-71, 2.2E-94, 4E-117, 7.4E-140}; // (same as loop above, rounded)
    }
    //----------

    flux_y = vector<vector<double> >(ngrd+1, vector<double>(nsp, 0.0));
    if(doSoot) flux_soot = vector<vector<double> >(ngrd+1, vector<double>(nsoot, 0.0));
    flux_h.resize(ngrd+1);

    //---------- radiation object

    doRadiation = false;

    radProps = make_shared<rad_planck_mean>();  // set even if doRadiation is false, since we switch it on/off for some cases

    Ttarget = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Set spatial grid (x and dx)
/// @param _L \input domain length (resets class data member L, so cases can be run with different L)
///
/// |  *  |  *  |  *  |  *  |  *  |
///   dx     dx    dx    dx    dx
///
////////////////////////////////////////////////////////////////////////////////

void fuego::setGrid(double _L) {

    L = _L;
    dx = vector<double>(ngrd, L/ngrd);
    x.resize(ngrd);

    if(isPremixed) {       //----------- uniform grid
        x[0] = dx[0]/2;
        for (size_t i=1; i<ngrd; i++)
            x[i] = x[i-1] + (dx[i-1]+dx[i])/2;
    }
    else {                  //--------- segmented grid

        double Lfrac = 0.5;         // first Lfrac fraction of the domain length
        double Gfrac = 0.5;         // gets this Gfrac fraction of the grid points
        int n1 = ngrd*Gfrac;
        double dx1 = L*Lfrac/n1;
        for(int i=0; i<n1; i++)
            dx[i] = dx1;

        int n2 = ngrd-n1;
        double dx2 = L*(1.0-Lfrac)/n2;
        for(int i=n1; i<ngrd; i++)
            dx[i] = dx2;

        x[0] = dx[0]/2;
        for (size_t i=1; i<ngrd; i++)
            x[i] = x[i-1] + (dx[i-1]+dx[i])/2;
    }

    //--------- set fx
    // interpolation factors to interior faces
    // 0     1                      2           3     4
    // |  *  |           *          |     *     |  *  |
    //    0              1                2        3
    // vf[1] = v[0]*(dx[1]/(dx[0]+dx[1])) + v[1]*(dx[0]/(dx[0]+dx[1]))
    //       = v[0]*(       fl[0]       ) + v[1]*(    1.0 - fl[0]    )
    //       = v[0]*(       fl[0]       ) + v[1]*(       fr[0]       )
    // 
    // vf[3] = v[2]*(dx[3]/(dx[2]+dx[3])) + v[3]*(dx[2]/(dx[2]+dx[3]))
    //       = v[2]*(       fl[2]       ) + v[3]*(    1.0 - fl[2]    )
    //       = v[2]*(       fl[2]       ) + v[3]*(       fr[2]       )
    //
    // vf[i] = v[i-1]*fl[i-1] + v[i]*fr[i-1]

    fl.resize(ngrd-1);
    fr.resize(ngrd-1);
    for(size_t i=0; i<ngrd-1; i++) {
        fl[i] = dx[i+1]/(dx[i]+dx[i+1]);
        fr[i] = 1.0-fl[i];
    }

}

////////////////////////////////////////////////////////////////////////////////
/// 
///  Write output file
///  @param fname \input write output to this file name (includes path)
/// 
////////////////////////////////////////////////////////////////////////////////

void fuego::writeFile(const string fname) {

    //-------------- compute auxiliary quantities

    vector<double> mixf;
    double mixfLbc;
    double mixfRbc;
    if(!isPremixed) {
        mixf.resize(ngrd);
        for(int i=0; i<ngrd; i++)
            mixf[i] = strm->getMixtureFraction(&y[i][0]);
        mixfLbc = strm->getMixtureFraction(&yLbc[0]);
        mixfRbc = strm->getMixtureFraction(&yRbc[0]);
    }

    vector<double> rho(ngrd);
    vector<double> h(ngrd);
    for(int i=0; i<ngrd; i++) {
        gas->setState_TPY(T[i], P, &y[i][0]);
        rho[i] = gas->density();
        h[i] = gas->enthalpy_mass();
    }
    gas->setState_TPY(TLbc, P, &yLbc[0]);
    double rhoLbc = gas->density();
    gas->setState_TPY(TRbc, P, &yRbc[0]);
    double rhoRbc = gas->density();

    //-------------- 

    ofstream ofile(fname.c_str());
    if(!ofile) {
        cerr << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
        exit(0);
    }

    ofile << "# L (m)  = " << L << endl;
    ofile << "# P (Pa) = " << P << endl;

    ofile << "#";
    int j=1;
    ofile << setw(15) << "00" << j++ << "_x";
    if(!isPremixed) ofile << setw(13) << "00" << j++ << "_mixf";
    ofile << setw(16) << "00" << j++ << "_T";
    ofile << setw(16) << "00" << j++ << "_h";
    ofile << setw(10) << "00" << j++ << "_density";
    for(int k=0; k<nsp; k++) {
        stringstream ss; ss << setfill('0') << setw(3) << j++ << "_" << gas->speciesName(k);
        ofile << setw(19) << ss.str();
    }
    if(doSoot) { // todo: make sure this is calling the correct variables // jansenpb 
        for(int u=0; u<nsoot; u++) {
            stringstream yy; yy << setfill('0') << setw(3) << j++ << "_" << "M"<<u;
            ofile << setw(19) << yy.str();
    }}

    ofile << scientific;
    ofile << setprecision(10);

    ofile << endl;
    ofile << setw(19) << 0;
    if(!isPremixed) ofile << setw(19) << mixfLbc;
    ofile << setw(19) << TLbc;
    ofile << setw(19) << hLbc;
    ofile << setw(19) << rhoLbc;
    for(int k=0; k<nsp; k++)
        ofile << setw(19) << yLbc[k];
    if(doSoot) {
        for(int k=0; k<nsoot; k++)
            ofile << setw(19) << 0;
    }

    for(int i=0; i<ngrd; i++) {
        ofile << endl;
        ofile << setw(19) << x[i];
        if(!isPremixed) ofile << setw(19) << mixf[i];
        ofile << setw(19) << T[i];
        ofile << setw(19) << h[i];
        ofile << setw(19) << rho[i];
        for(int k=0; k<nsp; k++)
            ofile << setw(19) << y[i][k];
        if(doSoot) { 
            for(int k=0; k<nsoot; k++)
                ofile << setw(19) << sootvars[i][k];         // todo: make sure the right vars are called
        }
    }

    if(isPremixed) {
        ofile << endl;
        ofile << setw(19) << L;
        ofile << setw(19) << T.back();
        ofile << setw(19) << h.back();
        ofile << setw(19) << rho.back();
        for(int k=0; k<nsp; k++)
            ofile << setw(19) << y.back()[k];
        if(doSoot) {
            for(int k=0; k<nsoot; k++)
                ofile << setw(19) << sootvars.back()[k];
        }
    }
    else {
        ofile << endl;
        ofile << setw(19) << L;
        ofile << setw(19) << mixfRbc;
        ofile << setw(19) << TRbc;
        ofile << setw(19) << hRbc;
        ofile << setw(19) << rhoRbc;
        for(int k=0; k<nsp; k++)
            ofile << setw(19) << yRbc[k];
        if(doSoot) {
            for(int k=0; k<nsoot; k++)
                ofile << setw(19) << 0.0;
        }
    }

    ofile.close();
}

////////////////////////////////////////////////////////////////////////////////
///
/// Store system state so that we can use this state later as an initial condition.
/// See setIC
///
////////////////////////////////////////////////////////////////////////////////

void fuego::storeState() {

    Pstore = P;
    ystore = y;
    Tstore = T;
    sootstore = sootvars;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Set initial condition for solution (or initial guess for steady solvers).
/// @param icType \input string indicating the type of the initial condition to use
/// @param fname \input read input from this file name (defaults to empty string)
///
////////////////////////////////////////////////////////////////////////////////

void fuego::setIC(const std::string icType, string fname) {

    if (icType == "linear") {
        gas->setState_TPY(TLbc, P, &yLbc[0]);
        double hLbc = gas->enthalpy_mass();
        gas->setState_TPY(TRbc, P, &yRbc[0]);
        double hRbc = gas->enthalpy_mass();
        for(int i=0; i<ngrd; i++) {
            for(int k=0; k<nsp; k++)
                y[i][k] = yLbc[k] + x[i]/L*(yRbc[k]-yLbc[k]);
            double h = hLbc + x[i]/L*(hRbc-hLbc);
            gas->setMassFractions(&y[i][0]);
            gas->setState_HP(h,P);
            T[i] = doEnergyEqn ? gas->temperature() : LI->interp(x[i]);
        }
    }

    //-------------------

    else if (icType == "equilibrium") {
        setIC("linear");
        for(int i=0; i<ngrd; i++) {
            gas->setState_TPY(T[i], P, &y[i][0]);
            gas->equilibrate("HP");
            gas->getMassFractions(&y[i][0]);
            T[i] = doEnergyEqn ? gas->temperature() : LI->interp(x[i]);   // redundant
        }
    }

    //-------------------

    else if (icType == "stored") {
        P        = Pstore;
        y        = ystore;
        T        = Tstore;
        sootvars = sootstore;
    }

    //-------------------

    else if (icType == "premixed") {
        if(doEnergyEqn) {
            gas->setMassFractions(&yLbc[0]);
            gas->setState_HP(hLbc,P);
            gas->equilibrate("HP");
            for(int i=0; i<ngrd; i++) {
                gas->getMassFractions(&y[i][0]);
                T[i] = gas->temperature();
            }
        }
        else {
            for(int i=0; i<ngrd; i++) {
                gas->setState_TPY(LI->interp(x[i]), P, &yLbc[0]);
                gas->equilibrate("TP");
                gas->getMassFractions(&y[i][0]);
                T[i] = gas->temperature();
            }
        }
    }

    //-------------------

    // else if (icType == "file") {{{{
    //     ifstream ifile(fname.c_str());
    //     if(!ifile) {
    //         cerr << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
    //         exit(0);
    //     }
    //     string s1;
    //     getline(ifile, s1);       // get header lines
    //     getline(ifile, s1);
    //     getline(ifile, s1);

    //     vector<double> XX;
    //     vector<double> HH;
    //     vector<vector<double> > YY;
    //     double d;
    //     while(!ifile.eof()) {
    //         ifile >> d; XX.push_back(d);   // x
    //         ifile >> d;                    // T
    //         ifile >> d; HH.push_back(d);   // h
    //         YY.push_back(vector<double>());
    //         for(int k=0; k<nsp; k++) {
    //             ifile >> d; YY.back().push_back(d);    
    //         }
    //     }

    // }}}}

    //-------------------

    // Soot is inialized to zero in the constructor

    storeState();

}

////////////////////////////////////////////////////////////////////////////////
///
/// Compute fluxes for all transported variables assuming unity Lewis numbers for energy and species.
/// Soot uses a thermophoretic flux.
/// Only called if doLe1 is false (the default).
///
////////////////////////////////////////////////////////////////////////////////

void fuego::setFluxesUnity() {

    //---------- cell center density and diffusivity

    vector<double> D(ngrd);
    vector<double> density(ngrd);
    vector<double> nu; if(doSoot) nu.resize(ngrd);
    vector<double> h(ngrd);
    for(int i=0; i<ngrd; i++) {
        gas->setMassFractions_NoNorm(&y[i][0]);
        gas->setState_TP(T[i], P);
        density[i] = gas->density();
        D[i] = trn->thermalConductivity()/(density[i]*gas->cp_mass());
        h[i] = gas->enthalpy_mass();
        if(doSoot) nu[i] = trn->viscosity()/density[i];
    }

    //---------- interpolate to face center density and diffusivity

    vector<double> D_f(ngrd+1);
    vector<double> density_f(ngrd+1);
    vector<double> nu_f; if(doSoot) nu_f.resize(ngrd+1);
    vector<double> T_f;  if(doSoot) T_f.resize(ngrd+1);

    gas->setState_TPY(TLbc, P, &yLbc[0]);       // this is only approximate for composition for premixed burner
    density_f[0] = gas->density();
    D_f[0] = trn->thermalConductivity()/(density_f[0]*gas->cp_mass());
    if(doSoot){
        nu_f[0] = trn->viscosity()/density_f[0];
        T_f[0] = TLbc;
    }

    if(!isPremixed) {
        gas->setState_TPY(TRbc, P, &yRbc[0]);
        density_f.back() = gas->density();
        D_f.back() = trn->thermalConductivity()/(density_f.back()*gas->cp_mass());
        if(doSoot) {
            nu_f.back() = trn->viscosity()/density_f.back();
            T_f.back() = TRbc;
        }
    }
    else{                          // extrapolate half a cell
        density_f.back() = density.back() +(density[ngrd-1] - density[ngrd-2])/(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
        D_f.back() = D.back() +(D[ngrd-1] - D[ngrd-2])/(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
        if(doSoot){ 
            nu_f.back() = nu.back() +(nu[ngrd-1] - nu[ngrd-2])/(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
            T_f.back() = T.back() +(T[ngrd-1] - T[ngrd-2])/(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
        }
    }

    for(int i=1, im=0; i<ngrd; i++, im++) {    // interpolate
        density_f[i] = density[im]*fl[im] + density[i]*fr[im];
        D_f[i]       = D[im]      *fl[im] + D[i]      *fr[im];
        if(doSoot) {
            nu_f[i] = nu[im]*fl[im] + nu[i]*fr[im]; 
            T_f[i]  = T[im] *fl[im] + T[i] *fr[im]; 
        }
    }

    //---------- fluxes y

    for(int k=0; k<nsp; k++) {
        if(isPremixed) {
            flux_y[0][k]    = mflux * yLbc[k];
            flux_y[ngrd][k] = mflux * y[ngrd-1][k];
        }
        else {
            flux_y[0][k]    = -density_f[0]    *D_f[0]    *(y[0][k]-yLbc[k])     *2/dx[0];
            flux_y[ngrd][k] = -density_f.back()*D_f.back()*(yRbc[k]-y[ngrd-1][k])*2/dx.back();
        }
        for(int i=1; i<ngrd; i++) {
            flux_y[i][k] = -density_f[i]*D_f[i]*(y[i][k]-y[i-1][k])*2/(dx[i-1]+dx[i]);
            if(isPremixed) flux_y[i][k] += mflux*y[i-1][k];
        }
    }

    //---------- fluxes h

    for(int i=1; i<ngrd; i++) {
        flux_h[i] = -density_f[i]*D_f[i]*(h[i]-h[i-1])*2/(dx[i-1]+dx[i]);
        if(isPremixed) flux_h[i] += mflux*h[i-1];
    }
    if(!isPremixed) {
        flux_h.back() = -density_f.back()*D_f.back()*(hRbc-h.back())*2/dx.back(); 
        flux_h[0]     = -density_f[0]*D_f[0]*(h[0]-hLbc)*2/dx[0];
    }
    else {
        flux_h.back() = mflux*h[ngrd-1];

        flux_h[0] = 0.0;
        vector<double> hsp(nsp);                     // J/kg species i
        gas->getEnthalpy_RT(&hsp[0]);                // (hhat/RT) where hhat = J/kmol
        for(size_t k=0; k<nsp; k++){                  // --> hsp = J/kg species i
            hsp[k] *= TLbc*Cantera::GasConstant/gas->molecularWeight(k);
            flux_h[0] += hsp[k]*flux_y[0][k];
        }
        gas->setState_TPY(TLbc, P, &yLbc[0]);
        flux_h[0] -= trn->thermalConductivity() * (T[1]-TLbc)/(dx[0]*0.5);
    }

    //---------- fluxes soot

    if(doSoot) {                // thermophoretic
        for(int k=0; k<nsoot; k++) {
            if(isPremixed) {
                flux_soot[0][k]    = 0.0;
                flux_soot[ngrd][k] = mflux/density_f.back()*sootvars[ngrd-1][k] - 
                                     0.556*sootvars[ngrd-1][k]*nu_f.back()/T_f.back()*
                                     (T_f.back()-T.back())/dx.back()*2;
            }
            else {
                flux_soot[0][k]    = -0.556*sootvars[0][k]*nu_f[0]/TLbc*(T[0]-TLbc)/dx[0]*2;
                flux_soot[ngrd][k] = -0.556*sootvars[ngrd-1][k]*nu_f.back()/TRbc*(TRbc-T.back())/dx.back()*2;
            }
            for(int i=1; i<ngrd; i++) {
                flux_soot[i][k] = 0.556*nu_f[i]*(T[i]-T[i-1])/T_f[i]*2/(dx[i-1]+dx[i]);
                flux_soot[i][k] *= (flux_soot[i][k] > 0 ? -sootvars[i][k] : -sootvars[i-1][k]); // upwind
                if(isPremixed) flux_soot[i][k] += mflux/density_f[i] * sootvars[i-1][k];        // upwind
            }
        }
    }

     /*if(doSoot) {            // unity Le like species
         for(int k=0; k<nsoot; k++) {
             if(isPremixed) {
                 flux_soot[0][k]    = 0.0;
                 flux_soot[ngrd][k] = mflux/density_f.back() * sootvars[ngrd-1][k];
             }
             else {
                 flux_soot[0][k]    = -density_f[0]    *D_f[0]    *(sootvars[0][k]/density[0]-0.0)     *2/dx[0];
                 flux_soot[ngrd][k] = -density_f.back()*D_f.back()*(0.0-sootvars[ngrd-1][k]/density[ngrd-1])*2/dx.back();
             }
             for(int i=1; i<ngrd; i++) {
                 flux_soot[i][k] = -density_f[i]*D_f[i]*(sootvars[i][k]/density[i]-sootvars[i-1][k]/density[i-1])*2/(dx[i-1]+dx[i]);
                 if(isPremixed) flux_soot[i][k] += mflux/density_f[i] * sootvars[i-1][k];
             }
         }
     }*/
}

////////////////////////////////////////////////////////////////////////////////
///
/// Compute fluxes for all transported variables assuming mixture-averaged transport for species and energy.
/// Soot uses a thermophoretic flux.
/// Only called if doLe1 is false (the default).
///
////////////////////////////////////////////////////////////////////////////////

void fuego::setFluxes() {

    //---------- cell center density and diffusivity

    vector<vector<double>> D(ngrd, vector<double> (nsp,0.0));
    vector<double> density(ngrd);
    vector<double> M(ngrd);
    vector<double> tcond(ngrd);
    vector<double> nu; if(doSoot) nu.resize(ngrd);
    for(int i=0; i<ngrd; i++) {
        gas->setMassFractions_NoNorm(&y[i][0]);
        gas->setState_TP(T[i], P);
        density[i] = gas->density();
        M[i] = gas->meanMolecularWeight();
        trn->getMixDiffCoeffs(&D[i][0]);
        tcond[i] = trn->thermalConductivity();
        if(doSoot)
            nu[i] = trn->viscosity()/density[i];
    }

    //---------- interpolate to face center density and diffusivity

    vector<vector<double>> D_f(ngrd+1, vector<double> (nsp,0.0));
    vector<vector<double>> y_f(ngrd+1, vector<double> (nsp,0.0));
    vector<double> density_f(ngrd+1);
    vector<double> M_f(ngrd+1);
    vector<double> T_f(ngrd+1);
    vector<double> tcond_f(ngrd+1);
    vector<double> nu_f; if(doSoot) nu_f.resize(ngrd+1);

    gas->setState_TPY(TLbc, P, &yLbc[0]);   // this is only approximate for composition for premixed burner
    density_f[0] = gas->density();
    M_f[0] = gas->meanMolecularWeight();
    trn->getMixDiffCoeffs(&D_f[0][0]);
    tcond_f[0] = trn->thermalConductivity();
    if(doSoot) nu_f[0]    = trn->viscosity()/density_f[0];
    T_f[0] = TLbc;
    y_f[0] = yLbc;

    if(!isPremixed) {
        gas->setState_TPY(TRbc, P, &yRbc[0]);
        density_f.back() = gas->density();
        M_f.back() = gas->meanMolecularWeight();
        trn->getMixDiffCoeffs(&(D_f.back()[0]));
        tcond_f.back() = trn->thermalConductivity();
        T_f.back() = TRbc;
        y_f.back() = yRbc;
        if(doSoot)
            nu_f.back() = trn->viscosity()/density_f.back();
    }
    else {

        density_f.back()  = density.back() +(density[ngrd-1] - density[ngrd-2])/(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
        tcond_f.back()    = tcond.back()   +(tcond[ngrd-1]   - tcond[ngrd-2])  /(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
        T_f.back()        = T.back()       +(T[ngrd-1]       - T[ngrd-2])      /(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
        M_f.back()        = M.back()       +(M[ngrd-1]       - M[ngrd-2])      /(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
        for(int k=0; k<nsp; k++) {
            D_f.back()[k] = D.back()[k]    +(D[ngrd-1][k]    - D[ngrd-2][k])   /(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
            y_f.back()[k] = y.back()[k]    +(y[ngrd-1][k]    - y[ngrd-2][k])   /(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
        }
        if(doSoot)
            nu_f.back()   = nu.back()      +(nu[ngrd-1]      - nu[ngrd-2])     /(dx[ngrd-1]+dx[ngrd-2])*2.0*(L-x[ngrd-1]);
    }

    for (int i=1, im=0; i<ngrd; i++, im++) {
        density_f[i]  = density[im]*fl[im] + density[i]*fr[im];
        T_f[i]        = T[im]      *fl[im] + T[i]      *fr[im];
        tcond_f[i]    = tcond[im]  *fl[im] + tcond[i]  *fr[im];
        M_f[i]        = M[im]      *fl[im] + M[i]      *fr[im];
        if(doSoot)
            nu_f[i]   = nu[im]     *fl[im] + nu[i]     *fr[im];
        for(int k=0; k<nsp; k++) {
            D_f[i][k] = D[im][k]   *fl[im] + D[i][k]   *fr[im];
            y_f[i][k] = y[im][k]   *fl[im] + y[i][k]   *fr[im];
        }
    }

    //---------- fluxes y

    double jstar;             // correction flux so that all fluxes sum to zero. This is equal to using a correction velocity
                              // j_i_corrected = j_i - Yi*jstar; jstar = sum(j_i).

    //--------- Boundary faces

    // Left boundary

    if(isPremixed)
        for(int k=0; k<nsp; k++)
            flux_y[0][k]    = mflux * yLbc[k];
    else {
        jstar = 0.0;
        for(int k=0; k<nsp; k++) {
            flux_y[0][k]    = -density_f[0]*D_f[0][k]*(y[0][k]-yLbc[k])*2/dx[0]
                              -density_f[0]*D_f[0][k]*y_f[0][k]*(M[0]-M_f[0])*2/dx[0]/M_f[0];
            jstar += flux_y[0][k];
        }	
        for(int k=0; k<nsp; k++) {
            flux_y[0][k] -= y_f[0][k]*jstar;
        }
    }

    // Right boundary

    if(isPremixed)
        for(int k=0; k<nsp; k++)
            flux_y[ngrd][k]    = mflux * y[ngrd-1][k];
    else {
        jstar = 0.0;
        for(int k=0; k<nsp; k++) {
            flux_y[ngrd][k] = -density_f.back()*D_f.back()[k]*(yRbc[k]-y[ngrd-1][k])*2/dx.back()
                              -density_f.back()*D_f.back()[k]*y_f.back()[k]*(M_f.back()-M[ngrd-1])*2/dx.back()/M_f.back();
            jstar += flux_y[ngrd][k];
        }
        for(int k=0; k<nsp; k++) {
            flux_y[ngrd][k] -= y_f[ngrd][k]*jstar;
        }
    }

    // Interior

    for (int i=1; i<ngrd; i++) {
        jstar = 0.0;
        for(int k=0; k<nsp; k++) {
            flux_y[i][k] = -density_f[i]*D_f[i][k]*(y[i][k]-y[i-1][k])*2/(dx[i-1]+dx[i])
                           -density_f[i]*D_f[i][k]*y_f[i][k]*(M[i]-M[i-1])*2/(dx[i-1]+dx[i])/M_f[i];
            jstar += flux_y[i][k];
        }
        for(int k=0; k<nsp; k++) {
            flux_y[i][k] -= y_f[i][k]*jstar;
            if(isPremixed) flux_y[i][k] += mflux*y[i-1][k];
        }
    }

    //---------- fluxes h

    // thermal conductivity portion

    flux_h[0]     = -tcond_f[0]*(T[0]-T_f[0])*2/dx[0];
    flux_h.back() = -tcond_f.back()*(T_f.back()-T.back())*2/dx.back(); 
    for(int i=1; i<ngrd; i++)
        flux_h[i] = -tcond_f[i]*(T[i]-T[i-1])*2/(dx[i-1]+dx[i]);

    // species portion (builds in advective flux if isPremixed)

    vector<double> hsp(nsp);

    gas->setState_TPY(TLbc, P, &yLbc[0]);
    gas->getEnthalpy_RT(&hsp[0]);
    for(size_t k=0; k<nsp; k++)
        flux_h[0] += flux_y[0][k]*hsp[k]*TLbc*Cantera::GasConstant/gas->molecularWeight(k);

    gas->setState_TPY(T_f.back(), P, &(y_f.back()[0]));
    gas->getEnthalpy_RT(&hsp[0]);
    for(size_t k=0; k<nsp; k++)
        flux_h.back() += flux_y.back()[k]*hsp[k]*T_f.back()*Cantera::GasConstant/gas->molecularWeight(k);

    for(size_t i=1; i<ngrd; i++) {
        gas->setState_TP(T_f[i], P);        // hsp depends on T but not on y
        gas->getEnthalpy_RT(&hsp[0]);
        for(size_t k=0; k<nsp; k++)
            flux_h[i] += flux_y[i][k]*hsp[k]*T_f[i]*Cantera::GasConstant/gas->molecularWeight(k);
    }

    //---------- fluxes soot

    if(doSoot) {                // thermophoretic
        for(int k=0; k<nsoot; k++) {
            if(isPremixed) {
                flux_soot[0][k]    = 0.0;
                flux_soot[ngrd][k] = mflux/density_f.back()*sootvars[ngrd-1][k] - 
                                     0.556*sootvars[ngrd-1][k]*nu_f.back()/T_f.back()*
                                     (T_f.back()-T.back())/dx.back()*2;
                for(int i=1; i<ngrd; i++) {
                    flux_soot[i][k]  = -0.556*nu_f[i]*(T[i]-T[i-1])/T_f[i]*2/(dx[i-1]+dx[i]) + mflux/density_f[i]; // soot velocity
                    flux_soot[i][k] *= (flux_soot[i][k] > 0 ? sootvars[i-1][k] : sootvars[i][k]);                  // upwind flux
                }
            }
            else {
                flux_soot[0][k]    = -0.556*sootvars[0][k]*nu_f[0]/TLbc*(T[0]-TLbc)/dx[0]*2;
                flux_soot[ngrd][k] = -0.556*sootvars[ngrd-1][k]*nu_f.back()/TRbc*(TRbc-T.back())/dx.back()*2;
                for(int i=1; i<ngrd; i++) {
                    flux_soot[i][k]  = -0.556*nu_f[i]*(T[i]-T[i-1])/T_f[i]*2/(dx[i-1]+dx[i]);     // soot velocity
                    flux_soot[i][k] *= (flux_soot[i][k] > 0 ? sootvars[i-1][k] : sootvars[i][k]); // upwind
                }

                //flux_soot[0][k] += -nu_f[0]*(sootvars[0][k]-0.0)/dx[0]*2;
                //flux_soot[ngrd][k] += -nu_f.back()*(0.0-sootvars[ngrd-1][k])/dx.back()*2;
                //for(int i=1; i<ngrd; i++)
                //    flux_soot[i][k] += -nu_f[i]*(sootvars[i][k]-sootvars[i-1][k])/(dx[i-1]+dx[i])*2;
            }
        }
    }

}

////////////////////////////////////////////////////////////////////////////////
///
/// Solve steady state problem. Uses Sundials Kinsol.
/// It is more robust to solve the unsteady problem to steady state.
/// Solves F(vars) = 0, where vars are the vector of variables at all grid points 
///   and F is the equation for each of them. See fuego::Func.
///
////////////////////////////////////////////////////////////////////////////////

void fuego::solveSS() {

    //---------- transfer variables into single array

    vector<double> vars(nvarA);

    for(size_t i=0; i<ngrd; i++) {
        for(size_t k=0; k<nsp; k++)
            vars[Ia(i,k)] = y[i][k];
        vars[Ia(i,nvar-1)] = T[i]/Tscale;      // dolh comment to remove h
    }
    
    //-------------- set vars0, F0 for homotopy

    vars0 = vars;
    F0.resize(nvarA);
    Func(&vars0[0], &F0[0]);

    //---------- setup solver

    vector<double> scales_v(nvarA, 1.0);
    vector<double> scales_f(nvarA, 1.0);
    vector<double> constraints(nvarA, 0.0);

    for(int i=0; i<nvarA; i++) {
        scales_v[i] = vars[i] < 1E-3 ? 1E-3 : vars[i];
        scales_f[i] = vars[i] < 1E-3 ? 1E-3 : vars[i];
    }


    double ftol = 5E-6; //5E-4;   // default: 5E-6
    double stol = 2E-11;  // default: 2E-11
    int mu = nvar*2-1;
    int ml = nvar*2-1;
    solver_kinsol s_kin(this, nvarA, scales_v, scales_f, constraints, mu, ml, ftol, stol);

    //---------- solve

    // s = 0.0;
    // s_kin.solve(Func_kinsol, vars);
    // s = 0.5;
    // s_kin.solve(Func_kinsol, vars);

    for(s=0.0; s<=1.0; s+=0.01) {
        cout << endl << "s = " << s; cout.flush();
        s_kin.solve(Func_kinsol, vars);
    }

    //---------- transfer variables back

    for(size_t i=0; i<ngrd; i++) {
        for(size_t k=0; k<nsp; k++)
            y[i][k] = vars[Ia(i,k)];
        T[i] = vars[Ia(i,nvar-1)]*Tscale;      // dolh comment to remove h
    }
}

////////////////////////////////////////////////////////////////////////////////
///
/// System of equations solved by the steady solver, F(vars) = 0
/// @param vars \input pointer to array of variables at every grid point.
/// @param F    \output pointer to array of function values for each variable at each grid point
/// 1D vectors are accessed by convenience function Ia(igrid, kvar), as in vars[Ia(i,k)].
///
////////////////////////////////////////////////////////////////////////////////

int fuego::Func(const double *vars, double *F) {

    //------------ transfer variables

    for(size_t i=0; i<ngrd; i++) {
        for(size_t k=0; k<nsp; k++)
            y[i][k] = vars[Ia(i,k)];
        T[i] = vars[Ia(i,nvar-1)]*Tscale;      // dolh comment to remove h
    }

    //------------ set function values

    if(doLe1)
        setFluxesUnity();
    else
        setFluxes();

    vector<double> Q(ngrd);
    if(doRadiation) setQrad(Q);

    vector<double> rr(nsp);

    for(size_t i=0; i<ngrd; i++) {
        gas->setMassFractions_NoNorm(&y[i][0]);
        gas->setState_TP(T[i], P);
        kin->getNetProductionRates(&rr[0]);
        for(size_t k=0; k<nsp; k++)
            F[Ia(i,k)]  = flux_y[i+1][k] - flux_y[i][k]
                         - rr[k]*gas->molecularWeight(k)*dx[i];
        if(doEnergyEqn) {
            F[Ia(i,nvar-1)] = (flux_h[i+1] - flux_h[i])/hscale; // dolh comment to remove h
            if(doRadiation) F[Ia(i,nvar-1)] -= Q[i]/hscale;
        }
        else
            F[Ia(i,nvar-1)] = 0.0; // dolh comment to remove h
    }

    //------------ augment for homotopy

    for(size_t i=0; i<nvarA; i++)
        F[i] = F[i] + (s-1.0)*F0[i];

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Kinsol interface function. 
/// Kinsol calls this function, which then calls user_data's Func.
/// @param varsKS    \input vector of all variables at all grid points (KS for KinSol).
/// @param fvec      \output vector of all variables at all grid points.
/// @param user_data \inout pointer to user data ("this" fuego object).
///
////////////////////////////////////////////////////////////////////////////////

int Func_kinsol(N_Vector varsKS, N_Vector fvec, void *user_data) {

    fuego *flm = static_cast<fuego *>(user_data);

    double *vars = N_VGetArrayPointer(varsKS);
    double *F = N_VGetArrayPointer(fvec);

    int rv = flm->Func(vars, F);

    return rv;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Set radiative source term for the energy equation.
/// @param Q \output source term at each grid point (W/m3)
///
////////////////////////////////////////////////////////////////////////////////

void fuego::setQrad(vector<double> &Q) {

    vector<double> kabs, awts;
    double fvsoot = 0.0;
    double xH2O, xCO2, xCO, xCH4;
    int isp;

    for(int i=0; i<ngrd; i++) {
        isp = gas->speciesIndex("H2O"); 
        xH2O = y[i][isp]/gas->molecularWeight(isp)*gas->meanMolecularWeight();
        isp = gas->speciesIndex("CO2"); 
        xCO2 = y[i][isp]/gas->molecularWeight(isp)*gas->meanMolecularWeight();
        isp = gas->speciesIndex("CO"); 
        xCO = y[i][isp]/gas->molecularWeight(isp)*gas->meanMolecularWeight();
        isp = gas->speciesIndex("CH4"); 
        xCH4 = y[i][isp]/gas->molecularWeight(isp)*gas->meanMolecularWeight();
        radProps->get_k_a(kabs, awts, T[i], P, fvsoot, xH2O, xCO2, xCO, xCH4);

        Q[i] = -4.0*rad::sigma*kabs[0]*(pow(T[i],4.0) - pow(TLbc,4.0));
    }
}

////////////////////////////////////////////////////////////////////////////////
///
/// Solve unsteady fuego problems.
/// Assumes y, T are initialized
/// Two modes: write on every time step of size dt, or write on temperature steps of size dT.
/// Default is in time --> Tmin, Tmax are zero --> dT = 0
/// 
/// @param nTauRun     \input number of characteristic times to solve for.
/// @param nSteps      \input number of steps to take during the solution.
/// @param doWriteTime \input if true (default) then write solution in time
/// @param Tmin        \input minimum temperature (default is 0, see code) (for stepping to desired T)
/// @param Tmax        \input maximum temperature (default is 0, see code) (for stepping to desired T)
///
////////////////////////////////////////////////////////////////////////////////

void fuego::solveUnsteady(const double nTauRun, const int nSteps, const bool doWriteTime, 
                          const double Tmin, const double Tmax) {

    //---------- transfer variables into single array
    vector<double> vars(nvarA);

    for(size_t i=0; i<ngrd; i++) {
        for(size_t k=0; k<nsp; k++)
            vars[Ia(i,k)] = y[i][k];
        if(doSoot) {
            for(size_t k=nsp; k<nsp+nsoot; k++) {  // jansenpb
                vars[Ia(i,k)] = sootvars[i][k-nsp]/sootScales[k-nsp];    // dont't forget sootscales
            }
        }
        vars[Ia(i,nvar-1)] = T[i]/Tscale;      // dolh comment to remove h
    }
    
    //---------- setup solver

    double rtol = 1E-4;                        // 1E-4 --> error at 0.01%, cvode documentation recomends < 1E-3
    vector<double> atol(nvarA, 1.0);           // noise level of the given variable

    for(int i=0; i<nvarA; i++)
        atol[i] = 1E-12;
        //atol[i] = vars[i] < 1E-3 ? 1E-3 : vars[i];


    int mu = nvar*2-1;
    int ml = nvar*2-1;
    integrator_cvode integ(rhsf_cvode, this, nvarA, rtol, atol, mu, ml, vars);

    //---------- solve

    double D = 0.00035;                                   // avg thermal diffusivity
    double t = 0.0;
    if(isPremixed) gas->setState_TPY(TLbc, P, &yLbc[0]); // needed fro density in tend

    double tau;
    if(isPremixed)
        tau = L/(mflux/gas->density());
    else if(isFlamelet)
        tau = 1.0/chi0;
    else
        tau = L*L/D;
    double tend = nTauRun * tau;

    double dt = tend/nSteps;
    dT = Tmax==Tmin ? 0.0 : (Tmax-(Tmin+0.1))/nSteps;
    Ttarget = Tmax - dT;
    isave = 1;

    for(int istep=1; istep<=nSteps; istep++, t+=dt) {
        integ.integrate(vars, dt);
        if(doWriteTime && dT <= 0.0) {           // write in time; (write in Temp is in rhsf)
            stringstream ss; 
            if(isFlamelet)
                ss << "X_" << chi0 << "U_" << setfill('0') << setw(3) << isave++ << ".dat";
            else
                ss << "L_" << L    << "U_" << setfill('0') << setw(3) << isave++ << ".dat";
            string fname = ss.str();
            writeFile(fname);
        }
    }

    //---------- transfer variables back

    for(size_t i=0; i<ngrd; i++) {
        for(size_t k=0; k<nsp; k++)
            y[i][k] = vars[Ia(i,k)];
        if(doSoot) {
            for(size_t k=nsp; k<nsp+nsoot; k++)
                sootvars[i][k-nsp] = vars[Ia(i,k)]*sootScales[k-nsp];    // jansenpb
        }
        T[i] = vars[Ia(i,nvar-1)]*Tscale;      // dolh comment to remove h
    }

    Ttarget = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Right hand side function for unsteady solver, as in dvar/dt = rhsf(var)
/// @param vars    \input current value of variables at each grid point
/// @param dvarsdt \output rate of each variable at each grid point
///
////////////////////////////////////////////////////////////////////////////////

int fuego::rhsf(const double *vars, double *dvarsdt) {

    //------------ transfer variables

    for(size_t i=0; i<ngrd; i++) {
        for(size_t k=0; k<nsp; k++)
            y[i][k] = vars[Ia(i,k)];
        if(doSoot) {
            for(size_t k=nsp; k<nsp+nsoot; k++) {
                sootvars[i][k-nsp] = vars[Ia(i,k)]*sootScales[k-nsp];    // jansenpb
            }
        }
        T[i] = vars[Ia(i,nvar-1)]*Tscale;      // dolh comment to remove h
    }

    //------------ set rates dvarsdt (dvars/dt)

    if(doLe1)
        setFluxesUnity();
    else
        setFluxes();

    vector<double> rr(nsp);
    vector<double> Q(ngrd);
    if(doRadiation) setQrad(Q);

    vector<double> yPAH; if(doSoot) yPAH.resize(6,0.0);

    vector<double> yGasForSM; if(doSoot) yGasForSM.resize( (size_t)gasSp::size );

    for(size_t i=0; i<ngrd; i++) {
        gas->setMassFractions_NoNorm(&y[i][0]);
        gas->setState_TP(T[i], P);
        kin->getNetProductionRates(&rr[0]);          // kmol/m3*s
        double rho = gas->density();
        double mu  = trn->viscosity();
        if(doSoot) {
            for(size_t kSootGases=0; kSootGases<(size_t)gasSp::size; kSootGases++) {
                size_t kgas = gas->speciesIndex(gasSpMapIS[kSootGases]);
                yGasForSM[kSootGases] = (kgas != Cantera::npos) ? y[i][kgas] : 0.0;
            }
            SMstate->setState(T[i], P, rho, mu, yGasForSM, yPAH, sootvars[i], nsoot);
            SM->setSourceTerms(*SMstate);
        }
        for(size_t k=0; k<nsp; k++)
            dvarsdt[Ia(i,k)]  = -(flux_y[i+1][k] - flux_y[i][k])/(rho*dx[i]) + 
                                rr[k]*gas->molecularWeight(k)/rho;
        if(doSoot) {
            for(size_t k=nsp; k<nsp+nsoot; k++) {
                dvarsdt[Ia(i,k)] = -(flux_soot[i+1][k-nsp] - flux_soot[i][k-nsp])/(dx[i]) + SM->sources.sootSources[k-nsp];
                dvarsdt[Ia(i,k)] /= sootScales[k-nsp];     //jansenpb; to match Tscale below; line 913 
            }
            // loop over the gas species in the soot model and compare with Cantera
            // update the gas source terms from the soot model

            for(size_t kSootGases=0; kSootGases<(size_t)gasSp::size; kSootGases++) {
                size_t kgas = gas->speciesIndex(gasSpMapIS[kSootGases]);
                if(kgas != Cantera::npos)
                    dvarsdt[Ia(i,kgas)] += SM->sources.gasSources[kSootGases];
            }
        }

        if(doEnergyEqn) {
            double cp  = gas->cp_mass();
            vector<double> hsp(nsp);                     // J/kg species i
            gas->getEnthalpy_RT(&hsp[0]);                // (hhat/RT) where hhat = J/kmol
            for(size_t k=0; k<nsp; k++)                  // --> hsp = J/kg species i
                hsp[k] *= T[i]*Cantera::GasConstant/gas->molecularWeight(k);
            double sum_hkdykdt = 0.0;
            for(size_t k=0; k<nsp; k++)
                sum_hkdykdt += hsp[k]*dvarsdt[Ia(i,k)];
            dvarsdt[Ia(i,nvar-1)]  =  -(flux_h[i+1] - flux_h[i])/(rho*cp*dx[i]) - sum_hkdykdt/cp;
            if(doRadiation) dvarsdt[Ia(i,nvar-1)] += Q[i]/(rho*cp);
            dvarsdt[Ia(i,nvar-1)] /= Tscale;
        }
        else
            dvarsdt[Ia(i,nvar-1)] = 0.0;
    }
    //-------------
    double TmaxLocal = *max_element(T.begin(), T.end());
    if(TmaxLocal <= Ttarget) {
        cout << endl << isave << "  " << TmaxLocal << "  " << Ttarget << "  ";
        stringstream ss; ss << "L_" << L << "U_" << setfill('0') << setw(3) << isave++ << ".dat";
        string fname = ss.str();
        writeFile(fname);
        Ttarget -= dT;
    }

    //-------------

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Right hand side function for unsteady solver for flamelet equations, as in dvar/dt = rhsf(var)
/// @param vars    \input current value of variables at each grid point
/// @param dvarsdt \output rate of each variable at each grid point
///
////////////////////////////////////////////////////////////////////////////////

int fuego::rhsf_flamelet(const double *vars, double *dvarsdt) {

    //------------ transfer variables

    for(size_t i=0; i<ngrd; i++) {
        for(size_t k=0; k<nsp; k++)
            y[i][k] = vars[Ia(i,k)];
        if(doSoot)
            for(size_t k=nsp; k<nsp+nsoot; k++)
                sootvars[i][k-nsp] = vars[Ia(i,k)]*sootScales[k-nsp];    // jansenpb
        T[i] = vars[Ia(i,nvar-1)]*Tscale;      // dolh comment to remove h
    }

    //------ uncomment to solve for h instead (1 of 2)
    // for(size_t i=0; i<ngrd; i++) {
    //     gas->setMassFractions(&y[i][0]);
    //     gas->setState_HP(hLbc*(1.0-x[i]) + hRbc*x[i], P);
    //     T[i] = gas->temperature();
    // }

    //------------ set variables

    vector<vector<double> > rr(ngrd, vector<double>(nsp));
    vector<double>          rho(ngrd);

    vector<double>          cp(ngrd);
    vector<vector<double> > hsp(ngrd, vector<double>(nsp));
    vector<double>          hsprrSum(ngrd, 0.0);

    for(size_t i=0; i<ngrd; i++) {

        //------- set gas

        gas->setMassFractions_NoNorm(&y[i][0]);
        gas->setState_TP(T[i], P);

        //------- reaction rates

        kin->getNetProductionRates(&rr[i][0]);          // kmol/m3*s
        for(size_t k=0; k<nsp; k++)
            rr[i][k] *= gas->molecularWeight(k);        // kg/m3*2

        //------- species enthalpies

        gas->getEnthalpy_RT(&hsp[i][0]);                // (h/RT) where h = J/kmol
        for(size_t k=0; k<nsp; k++) {                   // --> hsp = J/kg species i
            hsp[i][k] *= T[i]*Cantera::GasConstant/gas->molecularWeight(k);
            hsprrSum[i] += hsp[i][k] * rr[i][k];
        }

        //-------- 

        rho[i] = gas->density();                        // kg/m3
        cp[i]  = gas->cp_mass();
    }

    //------------ species

    vector<double> d2ydz2(ngrd);
    vector<double> yy(ngrd);                            // intermediate transfer array
    for(size_t k=0; k<nsp; k++) {
        for(size_t i=0; i<ngrd; i++) 
            yy[i] = y[i][k];
        setDerivative2(yLbc[k], yRbc[k], yy, d2ydz2);
        for(size_t i=0; i<ngrd; i++)
            dvarsdt[Ia(i,k)] = 0.5*chi[i]*d2ydz2[i] + rr[i][k]/rho[i];
    }

    //------------ energy (temperature)

    vector<double> d2Tdz2(ngrd);
    setDerivative2(TLbc, TRbc, T, d2Tdz2);

    vector<double> dTdz(ngrd);
    setDerivative(TLbc, TRbc, T, dTdz);

    vector<double> dcpdz(ngrd);
    setDerivative(cpLbc, cpRbc, cp, dcpdz);

    vector<double> dydzdhdzSum(ngrd, 0.0);    //todo: fill this in
    vector<double> dykdz(ngrd);
    vector<double> dhkdz(ngrd);
    vector<double> hh(ngrd);                            // intermediate transfer array
    for(size_t k=0; k<nsp; k++) {
        for(size_t i=0; i<ngrd; i++) {
            yy[i] = y[i][k];
            hh[i] = hsp[i][k];
        }
        setDerivative(yLbc[k],   yRbc[k],   yy, dykdz);
        setDerivative(hspLbc[k], hspRbc[k], hh, dhkdz);
        for(size_t i=0; i<ngrd; i++)
            dydzdhdzSum[i] += dykdz[i]*dhkdz[i];
    }

    for(size_t i=0; i<ngrd; i++) {
        dvarsdt[Ia(i,nvar-1)] = -hsprrSum[i]/(cp[i]*rho[i]) + 0.5*chi[i]*
                                (d2Tdz2[i] + (dTdz[i]*dcpdz[i] + dydzdhdzSum[i])/cp[i]);
        dvarsdt[Ia(i,nvar-1)] /= Tscale;
    }

    //------ uncomment to solve for h instead (2 of 2)
    // for(size_t i=0; i<ngrd; i++)
    //     dvarsdt[Ia(i,nvar-1)] = 0.0;

    //---------- soot


    //-------------

    double TmaxLocal = *max_element(T.begin(), T.end());
    if(TmaxLocal <= Ttarget) {
        cout << endl << isave << "  " << TmaxLocal << "  " << Ttarget << "  ";
        stringstream ss; ss << "X_" << chi0 << "U_" << setfill('0') << setw(3) << isave++ << ".dat";
        string fname = ss.str();
        writeFile(fname);
        Ttarget -= dT;
    }

    //-------------

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Compute derivative of a profile at cell centers
/// @param vL   \input variable on left face of domain
/// @param vR   \input variable on right face of domain
/// @param v    \input variable grid points
/// @param dvdx \output derivative on grid points
///
/// Interpolate variable to faces, then take the derivative between faces
///       |  *  |           *          |     *     |  *  |
///
////////////////////////////////////////////////////////////////////////////////

void fuego::setDerivative(const double vL, const double vR, 
                          const vector<double> &v, vector<double> &dvdx) {

    double vfL = vL;                            // first cell (left side)
    double vfR = v[0]*fl[0]+v[1]*fr[0];
    dvdx[0] = (vfR-vfL)/dx[0];

    for(size_t i=1; i<ngrd-1; i++) {            // interior cells
        vfL = vfR;
        vfR = v[i]*fl[i]+v[i+1]*fr[i];
        dvdx[i] = (vfR-vfL)/dx[i];
    }
    vfL = vfR;                                  // last cell (right side)
    vfR = vR;
    dvdx[ngrd-1] = (vfR-vfL)/dx[ngrd-1];
}

////////////////////////////////////////////////////////////////////////////////
///
/// Compute second derivative of a profile at cell centers
/// @param vL   \input variable on left face of domain
/// @param vR   \input variable on right face of domain
/// @param v    \input variable grid points
/// @param d2vdx2 \output derivative on grid points
///
/// Get derivative on faces, then take difference between faces
///       |  *  |           *          |     *     |  *  |
///
////////////////////////////////////////////////////////////////////////////////

void fuego::setDerivative2(const double vL, const double vR, 
                           const vector<double> &v, vector<double> &d2vdx2) {

    double dvdxL = (v[0]-vL)  /dx[0]*2;               // first cell (left side)
    double dvdxR = (v[1]-v[0])/(dx[0]+dx[1])*2;
    d2vdx2[0]    = (dvdxR - dvdxL)/dx[0];

    for(size_t i=1; i<ngrd-1; i++) {                    // interior cells
        dvdxL     = dvdxR;
        dvdxR     = (v[i+1]-v[i])/(dx[i]+dx[i+1])*2;
        d2vdx2[i] = (dvdxR - dvdxL)/dx[i];
    }
    dvdxL = dvdxR;                                    // last cell (right side)
    dvdxR = (vR-v[ngrd-1])/dx[ngrd-1]*2;
    d2vdx2[ngrd-1] = (dvdxR - dvdxL)/dx[ngrd-1];
}

////////////////////////////////////////////////////////////////////////////////
///
/// Set scalar dissipation rate profile
/// @param _chi0 \input scalar dissipation rate multiplier
///
////////////////////////////////////////////////////////////////////////////////

void fuego::setChi(const double _chi0) {

    chi0 = _chi0;
    chi.resize(ngrd);
    double d, e;
    for(size_t i=0; i<ngrd; i++) {
        d = 2*x[i]-1;
        e = 1.0-d*d;
        chi[i] = chi0*e*e;
    }
}

////////////////////////////////////////////////////////////////////////////////
///
/// CVODE interface; CVODE calls this function, which then calls user_data's rhsf 
/// @param t         \input current time (not used here as there are no explicit time dependencies, like S(t)
/// @param varsCV    \input cvode variables (all vars at all grid points)
/// @param dvarsdtCV \output cvode rates of all variables (all vars at all grid points)
/// @param user_data \inout pointer to user data ("this" fuego object).
///
////////////////////////////////////////////////////////////////////////////////

int rhsf_cvode(realtype t, N_Vector varsCV, N_Vector dvarsdtCV, void *user_data) {
    fuego *flm = static_cast<fuego *>(user_data);

    double *vars  = N_VGetArrayPointer(varsCV);
    double *dvarsdt = N_VGetArrayPointer(dvarsdtCV);

    int rv;
    if(flm->isFlamelet)
        rv = flm->rhsf_flamelet(vars, dvarsdt);
    else
        rv = flm->rhsf(vars, dvarsdt);

    return rv;
}
