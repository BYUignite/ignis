#pragma once

#include <vector>

#include <cvode/cvode.h>               // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>    // access to serial N_Vector
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_band.h>  // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_band.h>  // access to dense SUNLinearSolver

////////////////////////////////////////////////////////////////////////////////
///
/// Inferface class for CVODE ODE integrator.
/// Header only
///
////////////////////////////////////////////////////////////////////////////////

class integrator_cvode {

public: 

//////////////////// DATA MEMBERS ////////////////////

    SUNContext      sun;           ///< sundials object
    void           *cmem;          ///< cvode object
    N_Vector        vars;          ///< vector of variables being solved
    N_Vector        atol;          ///< vector atol (absolute tolerance, for each variable)
    sunrealtype     rtol;          ///< scalar rtol (relative tolerance)
    unsigned        nvar;          ///< number of equations being solved
    SUNMatrix       J;             ///< matrix for linear solver
    SUNLinearSolver LS;            ///< linear solver
    int             rv;            ///< return value: checking status of calls

//////////////////// MEMBER FUNCTIONS ////////////////////

//////////////////////////////////////////////////////////
///
/// Constructor function
/// @param Func       \input pointer to right-hand-side function (rates) for the problem being solved
/// @param _user_data \inout point to user's class object with data and functions needed to compute the rates
/// @param _nvar      \input number of variables being solved
/// @param _rtol      \input relative tolerance for all variables
/// @param _atol      \input vector of absolute tolerances for all variables
/// @param mu         \input location of upper diagonal in the banded matrix for solution
/// @param ml         \input location of lower diagonal in the banded matrix for solution
/// @param y          \input initial condition of variables
///
//////////////////////////////////////////////////////////

integrator_cvode(
                 int (*Func)(sunrealtype, N_Vector, N_Vector, void*),
                 void *                _user_data,
                 const int             _nvar, 
                 const double          _rtol,
                 const std::vector<double> &_atol,
                 const int             mu,
                 const int             ml,
                 std::vector<double>  &y) :
                     nvar(_nvar),
                     rtol(_rtol) {

    rv = SUNContext_Create(0, &sun);

    vars = N_VNew_Serial(nvar, sun);
    atol = N_VNew_Serial(nvar, sun);
    for(int k=0; k<nvar; ++k) {
        NV_Ith_S(vars, k) = y[k];
        NV_Ith_S(atol, k) = _atol[k];
    }

    cmem = CVodeCreate(CV_BDF, sun);
    rv   = CVodeSetUserData(cmem, _user_data);
    rv   = CVodeInit(cmem, Func, 0.0, vars);
    rv   = CVodeSVtolerances(cmem, rtol, atol);
    rv   = CVodeSetMaxNumSteps(cmem, 20000); //5000);

    J    = SUNBandMatrix(nvar, mu, ml, sun);   // linear solver matrix J
    LS   = SUNLinSol_Band(vars, J, sun);          // set linear solver
    rv   = CVodeSetLinearSolver(cmem, LS, J);  // associate matrix J and solver LS
}

//////////////////////////////////////////////////////////
///
/// Main interface: integrate the ODE system
/// @param y  \inout vector of variables being solved (and initial condition)
/// @param dt \input time to intigrate for.
///
//////////////////////////////////////////////////////////

int integrate(std::vector<double> &y, const sunrealtype dt) {

    for(int k=0; k<nvar; ++k)
        NV_Ith_S(vars, k) = y[k];

    sunrealtype t;

    rv = CVodeReInit(cmem, 0.0, vars);
    rv = CVode(cmem, dt, vars, &t, CV_NORMAL);

    for(int k=0; k<nvar; ++k)
        y[k] = NV_Ith_S(vars,k);

    return rv;
}

//////////////////////////////////////////////////////////
///
/// Destructor, cleans up CVODE objects
///
//////////////////////////////////////////////////////////

~integrator_cvode() {

    N_VDestroy(vars);
    N_VDestroy(atol);
    CVodeFree(&cmem);
    SUNLinSolFree(LS);
    SUNMatDestroy(J);
    SUNContext_Free(&sun);
}
};
