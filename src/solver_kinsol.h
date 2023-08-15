#pragma once

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sunmatrix/sunmatrix_band.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_band.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */

#include <vector>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////

class solver_kinsol {


//////////////////// DATA MEMBERS //////////////////////

public: 

    void           *user_data;
    size_t          nvar;

    SUNContext      sun;
    void           *kmem;
    N_Vector        vars;
    N_Vector        scales_v;
    N_Vector        scales_f;
    N_Vector        constraints;
    SUNMatrix       J;
    SUNLinearSolver LS;
    int             rv;
    int             solver_type;                // KIN_NONE, KIN_LINESEARCH, KIN_FP, KIN_PICARD
    int             exact_or_modified_newton;   // 1 or 0, respectively

//////////////////// MEMBER FUNCTIONS /////////////////

    solver_kinsol(
        void                      *_user_data, 
        const size_t               _nvar,
        const std::vector<double> &_scales_v,
        const std::vector<double> &_scales_f,
        const std::vector<double> &_constraints,
        const int                  mu,
        const int                  ml,
        const double               ftol = 1E-5,
        const double               stol = 1E-5) :
            user_data(_user_data),
            nvar(_nvar) {

        solver_type              = KIN_LINESEARCH;
        exact_or_modified_newton = 1;

        rv = SUNContext_Create(NULL, &sun);

        vars        = N_VNew_Serial(nvar, sun);
        scales_v    = N_VNew_Serial(nvar, sun);
        scales_f    = N_VNew_Serial(nvar, sun);
        constraints = N_VNew_Serial(nvar, sun);

        // J  = SUNDenseMatrix(nvar, nvar, sun);   // linear solver matrix J
        J  = SUNBandMatrix(nvar, mu, ml,  sun);   // linear solver matrix J

        for(size_t k=0; k<nvar; k++) {
            NV_Ith_S(scales_v, k)    = _scales_v[k];
            NV_Ith_S(scales_f, k)    = _scales_f[k];
            NV_Ith_S(constraints, k) = _constraints[k];
        }

        kmem = KINCreate(sun);

        rv = KINSetUserData(     kmem, user_data);
        rv = KINSetConstraints(  kmem, constraints);
        rv = KINSetScaledStepTol(kmem, stol);
        rv = KINSetFuncNormTol(  kmem, ftol);
        rv = KINSetMaxSetupCalls(kmem, exact_or_modified_newton);
        rv = KINSetNumMaxIters(  kmem, 5000);
    }

//--------------------

int solve(int (*Func)(N_Vector, N_Vector, void*), 
          std::vector<double> &y) {

    //---------

    for(size_t k=0; k<nvar; k++)
        NV_Ith_S(vars, k) = y[k];

    rv = KINInit(kmem, Func, vars);
    // LS = SUNLinSol_Dense(vars, J, sun);           // set linear solver
    LS = SUNLinSol_Band(vars, J, sun);           // set linear solver
    rv = KINSetLinearSolver(kmem, LS, J);        // associate matrix J with solver LS

    //---------

    rv = KINSol(kmem, vars, solver_type, scales_v, scales_f);

    //---------

    for(size_t k=0; k<nvar; k++)
        y[k] = NV_Ith_S(vars,k);

    return rv;
}
//--------------

    ~solver_kinsol(){
        N_VDestroy(vars);
        N_VDestroy(scales_v);
        N_VDestroy(scales_f);
        N_VDestroy(constraints);
        KINFree(&kmem);
        SUNMatDestroy(J);
        SUNLinSolFree(LS);
        SUNContext_Free(&sun);
    }

};
