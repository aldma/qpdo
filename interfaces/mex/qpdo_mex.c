#include "mex.h"
#include <string.h>
#include "global_opts.h"
#include "qpdo.h"
#include "util.h"
#include "constants.h"
#include "cholmod.h"
#include "cholmod_matlab.h"
#include "cholmod_function.h"

//Modes of operation
#define MODE_DEFAULT_SETTINGS "default_settings"
#define MODE_RETRIEVE_CONSTANT "constant"
#define MODE_SETUP "setup"
#define MODE_WARM_START "warm_start"
#define MODE_UPDATE_BOUNDS "update_bounds"
#define MODE_UPDATE_LINEAR "update_q"
#define MODE_SOLVE "solve"
#define MODE_DELETE "delete"

// QPDO work identifier
static QPDOWorkspace* qpdo_work = NULL;

// all of the QPDO_INFO fieldnames as strings
const char* QPDO_INFO_FIELDS[] = {"iterations",     //c_int
                                  "oterations",     //c_int
                                  "status" ,        //char*
                                  "status_val" ,    //c_int
                                  "res_prim_norm",   //c_float
                                  "res_dual_norm",   //c_float
                                  "res_prim_in_norm",//c_float
                                  "res_dual_in_norm",//c_float
                                  "objective",      //c_float
                                  "setup_time",     //c_float, only used if PROFILING
                                  "solve_time",     //c_float, only used if PROFILING
                                  "run_time"};      //c_float, only used if PROFILING


const char* QPDO_SETTINGS_FIELDS[] = {"max_iter",                   //c_int
                                      "inner_max_iter",             //c_int
                                      "max_time",                   //c_float
                                      "eps_abs",                    //c_float
                                      "eps_abs_in",                 //c_float
                                      "eps_prim_inf",               //c_float
                                      "eps_dual_inf",               //c_float
                                      "rho",                        //c_float
                                      "theta",                      //c_float
                                      "delta",                      //c_float
                                      "mu_min",                     //c_float
                                      "proximal",                   //c_int
                                      "sigma_init",                 //c_float
                                      "sigma_upd",                  //c_float
                                      "sigma_min",                  //c_float
                                      "scaling",                    //c_int
                                      "reset_newton_iter",          //c_int
                                      "print_interval",             //c_int
                                      "verbose"};                   //c_int


// internal utility functions
void      castToDoubleArr(c_float *arr, double* arr_out, size_t len);
void      setToNaN(double* arr_out, size_t len);
void      copyMxStructToSettings(const mxArray*, QPDOSettings*);
mxArray*  copySettingsToMxStruct(QPDOSettings* settings);
mxArray*  copyInfoToMxStruct(QPDOInfo* info);

// Function that mex calls when it closes unexpectedly. Frees the workspace.
void exitFcn() {
  if (qpdo_work != NULL) {
      qpdo_cleanup(qpdo_work);
      qpdo_work = NULL;
  }
}

/**
 * The gateway function to QPDO
 *
 * qpdo_mex('default_settings');
 * qpdo_mex('constant', constant_name)
 * qpdo_mex('setup',n,m,Q,q,A,l,u,theSettings);
 * qpdo_mex('warm_start', x, y);
 * [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = qpdo_mex('solve');
 * qpdo_mex('delete');
 *
 * @param nlhs Number of output arguments
 * @param plhs Array of output argument pointers
 * @param nrhs Number of input arguments
 * @param prhs Array of input argument pointers
 */
void mexFunction(int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs []) {
    // Set function to call when mex closes unexpectedly
    mexAtExit(exitFcn);
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    // report the default settings
    if (strcmp(cmd, MODE_DEFAULT_SETTINGS) == 0) {
        // Warn if other commands were ignored
        if (nrhs > 2)
            mexWarnMsgTxt("Default settings: unexpected number of arguments.");
        // Create a Settings structure in default form and report the results
        // Useful for external solver packages (e.g. Yalmip) that want to
        // know which solver settings are supported
        QPDOSettings* defaults = (QPDOSettings *)mxCalloc(1,sizeof(QPDOSettings));
        qpdo_set_default_settings(defaults);
        plhs[0] = copySettingsToMxStruct(defaults);
        mxFree(defaults);
        return;
    } else if (strcmp(cmd, MODE_DELETE) == 0) {
        // clean up the problem workspace
        if(qpdo_work != NULL){
            qpdo_cleanup(qpdo_work);
            qpdo_work = NULL;
        }
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 1)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    } else if (strcmp(cmd, MODE_SETUP) == 0) {
        // throw an error if this is called more than once
        if(qpdo_work != NULL){
            mexErrMsgTxt("Solver is already initialized with problem data.");
        }
        // Create data and settings containers
        QPDOSettings* settings = (QPDOSettings *)mxCalloc(1,sizeof(QPDOSettings));
        QPDOData*     data     = (QPDOData *)mxCalloc(1,sizeof(QPDOData));
        // handle the problem data first.  Matlab-side
        // class wrapper is responsible for ensuring that
        // Q and A are sparse matrices,  everything
        // else is a dense vector and all inputs are
        // of compatible dimension
        // Get mxArrays
        const mxArray* Q  = prhs[3];
        const mxArray* q  = prhs[4];
        const mxArray* A  = prhs[5];
        const mxArray* l = prhs[6];
        const mxArray* u = prhs[7];
        // Create Data Structure
        data->n = (size_t) mxGetScalar(prhs[1]);
        data->m = (size_t) mxGetScalar(prhs[2]);
        data->c = 0; // TODO: Allow for a constant to be passed.
        data->q = mxGetPr(q);
        data->l = mxGetPr(l);
        data->u = mxGetPr(u);
        // Convert matrices from matlab to cholmod_sparse
        double dummy = 0;
        cholmod_sparse Amatrix, Qmatrix;
        data->A = sputil_get_sparse(A, &Amatrix, &dummy, 0);
        data->Q = sputil_get_sparse(Q, &Qmatrix, &dummy, -1);//Q is symmetric, use only lower part
        // Create Settings
        const mxArray* mxSettings = prhs[8];
        if(mxIsEmpty(mxSettings)){
            // use defaults
            qpdo_set_default_settings(settings);
        } else {
            // populate settings structure from mxArray input
            copyMxStructToSettings(mxSettings, settings);
        }
        // Setup workspace
        qpdo_work = qpdo_setup(data, settings);
        if(qpdo_work == NULL){
            mexErrMsgTxt("Invalid problem setup");
        }
         //cleanup temporary structures
         // Don't free data->q, data->l, data->l because they are pointers to the mxArrays
         // Don't free data->A and data->Q because they are only a shallow copy
         mxFree(data);
         mxFree(settings);
        return;
    } else if (strcmp(cmd, MODE_WARM_START) == 0) {
        if (nlhs != 0 || nrhs != 3 ){
            mexErrMsgTxt("Solve : wrong number of inputs / outputs");
        }
        if(!qpdo_work){
            mexErrMsgTxt("Work is not setup.");
        }
        c_float *x, *y;
        const mxArray* xmatlab  = prhs[1];
        const mxArray* ymatlab  = prhs[2];
        if (mxIsEmpty(xmatlab)) {
            x = NULL;
        } else {
            x = mxGetPr(xmatlab);
        }
        if (mxIsEmpty(ymatlab)) {
            y = NULL;
        } else {
            y = mxGetPr(ymatlab);
        }
        qpdo_warm_start(qpdo_work, x, y);
        return;
    } else if (strcmp(cmd, MODE_UPDATE_BOUNDS) == 0) {
        if (nlhs != 0 || nrhs != 3){
            mexErrMsgTxt("Update bounds : wrong number of inputs / outputs");
        }
        if(!qpdo_work){
            mexErrMsgTxt("Work is not setup.");
        }
        c_float *l, *u;
        const mxArray* l_matlab  = prhs[1];
        const mxArray* u_matlab  = prhs[2];
        if (mxIsEmpty(l_matlab)) {
            l = NULL;
        } else {
            l = mxGetPr(l_matlab);
        }
        if (mxIsEmpty(u_matlab)) {
            u = NULL;
        } else {
            u = mxGetPr(u_matlab);
        }
        qpdo_update_bounds(qpdo_work, l, u);
    } else if (strcmp(cmd, MODE_UPDATE_LINEAR) == 0) {
        if (nlhs != 0 || nrhs != 2){
            mexErrMsgTxt("Update q : wrong number of inputs / outputs");
        }
        if(!qpdo_work){
            mexErrMsgTxt("Work is not setup.");
        }
        if (!mxIsEmpty(prhs[1])) {
            c_float *q = mxGetPr(prhs[1]);
            qpdo_update_q(qpdo_work, q);
        } else {
            mexWarnMsgTxt("Update q: Empty q has no effect.");
        }
    } else if (strcmp(cmd, MODE_SOLVE) == 0) { // SOLVE
        if (nlhs != 5 || nrhs != 1){
            mexErrMsgTxt("Solve : wrong number of inputs / outputs");
        }
        if(!qpdo_work){
            mexErrMsgTxt("Work is not setup.");
        }
        qpdo_solve(qpdo_work);
        // Allocate space for solution
        // primal variables
        plhs[0] = mxCreateDoubleMatrix(qpdo_work->data->n,1,mxREAL);
        // dual variables
        plhs[1] = mxCreateDoubleMatrix(qpdo_work->data->m,1,mxREAL);
        // primal infeasibility certificate
        plhs[2] = mxCreateDoubleMatrix(qpdo_work->data->m,1,mxREAL);
        // dual infeasibility certificate
        plhs[3] = mxCreateDoubleMatrix(qpdo_work->data->n,1,mxREAL);
        // copy results to mxArray outputs
        // assume that five outputs will always
        // be returned to matlab-side class wrapper
        if ((qpdo_work->info->status_val != QPDO_PRIMAL_INFEASIBLE) &&
            (qpdo_work->info->status_val != QPDO_DUAL_INFEASIBLE)){
            // primal and dual solutions
            castToDoubleArr(qpdo_work->solution->x, mxGetPr(plhs[0]), qpdo_work->data->n);
            castToDoubleArr(qpdo_work->solution->y, mxGetPr(plhs[1]), qpdo_work->data->m);
            // infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), qpdo_work->data->m);
            setToNaN(mxGetPr(plhs[3]), qpdo_work->data->n);
        } else if (qpdo_work->info->status_val == QPDO_PRIMAL_INFEASIBLE){ //primal infeasible
            // primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), qpdo_work->data->n);
            setToNaN(mxGetPr(plhs[1]), qpdo_work->data->m);
            // primal infeasibility certificates
            castToDoubleArr(qpdo_work->dy, mxGetPr(plhs[2]), qpdo_work->data->m);
            // dual infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[3]), qpdo_work->data->n);
        } else if (qpdo_work->info->status_val == QPDO_DUAL_INFEASIBLE) { //dual infeasible
            // primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), qpdo_work->data->n);
            setToNaN(mxGetPr(plhs[1]), qpdo_work->data->m);
            // primal infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), qpdo_work->data->m);
            // dual infeasibility certificates
            castToDoubleArr(qpdo_work->dx, mxGetPr(plhs[3]), qpdo_work->data->n);
        } else if (qpdo_work->info->status_val == QPDO_PRIMAL_DUAL_INFEASIBLE){ // primal-dual infeasible
            // primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), qpdo_work->data->n);
            setToNaN(mxGetPr(plhs[1]), qpdo_work->data->m);
            // primal infeasibility certificates
            castToDoubleArr(qpdo_work->dy, mxGetPr(plhs[2]), qpdo_work->data->m);
            // dual infeasibility certificates
            castToDoubleArr(qpdo_work->dx, mxGetPr(plhs[3]), qpdo_work->data->n);
        }
        plhs[4] = copyInfoToMxStruct(qpdo_work->info); // Info structure
        return;
    } else if (strcmp(cmd, MODE_RETRIEVE_CONSTANT) == 0) { // Return solver constants
        char constant[32];
        mxGetString(prhs[1], constant, sizeof(constant));

        if (!strcmp("QPDO_INFTY", constant)){
            plhs[0] = mxCreateDoubleScalar(QPDO_INFTY);
            return;
        }
        if (!strcmp("QPDO_NAN", constant)){
            plhs[0] = mxCreateDoubleScalar(mxGetNaN());
            return;
        }
        if (!strcmp("QPDO_SOLVED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPDO_SOLVED);
            return;
        }
        if (!strcmp("QPDO_UNSOLVED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPDO_UNSOLVED);
            return;
        }
        if (!strcmp("QPDO_PRIMAL_INFEASIBLE", constant)){
            plhs[0] = mxCreateDoubleScalar(QPDO_PRIMAL_INFEASIBLE);
            return;
        }
        if (!strcmp("QPDO_DUAL_INFEASIBLE", constant)){
            plhs[0] = mxCreateDoubleScalar(QPDO_DUAL_INFEASIBLE);
            return;
        }
        if (!strcmp("QPDO_MAX_ITER_REACHED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPDO_MAX_ITER_REACHED);
            return;
        }
        mexErrMsgTxt("Constant not recognized.");
        return;
    } else {
        mexErrMsgTxt("Invalid QPDO mode");
    }
}

void castToDoubleArr(c_float *arr, double* arr_out, size_t len) {
    for (size_t i = 0; i < len; i++) {
        arr_out[i] = (double)arr[i];
    }
}

void setToNaN(double* arr_out, size_t len){
    size_t i;
    for (i = 0; i < len; i++) {
        arr_out[i] = mxGetNaN();
    }
}

mxArray* copyInfoToMxStruct(QPDOInfo* info){

    // create mxArray with the right number of fields
    int nfields  = sizeof(QPDO_INFO_FIELDS) / sizeof(QPDO_INFO_FIELDS[0]);
    mxArray* mxPtr = mxCreateStructMatrix(1,1,nfields,QPDO_INFO_FIELDS);

    // map the QPDO_INFO fields one at a time into mxArrays
    // matlab all numeric values as doubles
    mxSetField(mxPtr, 0, "iterations",        mxCreateDoubleScalar(info->iterations));
    mxSetField(mxPtr, 0, "oterations",        mxCreateDoubleScalar(info->oterations));
    mxSetField(mxPtr, 0, "status",            mxCreateString(info->status));
    mxSetField(mxPtr, 0, "status_val",        mxCreateDoubleScalar(info->status_val));
    mxSetField(mxPtr, 0, "res_prim_norm",      mxCreateDoubleScalar(info->res_prim_norm));
    mxSetField(mxPtr, 0, "res_dual_norm",      mxCreateDoubleScalar(info->res_dual_norm));
    mxSetField(mxPtr, 0, "res_prim_in_norm",   mxCreateDoubleScalar(info->res_prim_in_norm));
    mxSetField(mxPtr, 0, "res_dual_in_norm",   mxCreateDoubleScalar(info->res_dual_in_norm));
    mxSetField(mxPtr, 0, "objective",         mxCreateDoubleScalar(info->objective));

    #ifdef PROFILING
    //if not profiling, these fields will be empty
    mxSetField(mxPtr, 0, "setup_time",  mxCreateDoubleScalar(info->setup_time));
    mxSetField(mxPtr, 0, "solve_time",  mxCreateDoubleScalar(info->solve_time));
    mxSetField(mxPtr, 0, "run_time",    mxCreateDoubleScalar(info->run_time));
    #endif

    return mxPtr;

}

mxArray* copySettingsToMxStruct(QPDOSettings* settings){

    int nfields  = sizeof(QPDO_SETTINGS_FIELDS) / sizeof(QPDO_SETTINGS_FIELDS[0]);
    mxArray* mxPtr = mxCreateStructMatrix(1,1,nfields,QPDO_SETTINGS_FIELDS);

    // map the QPDO_SETTINGS fields one at a time into mxArrays
    // matlab handles everything as a double
    mxSetField(mxPtr, 0, "max_iter",                  mxCreateDoubleScalar(settings->max_iter));
    mxSetField(mxPtr, 0, "inner_max_iter",            mxCreateDoubleScalar(settings->inner_max_iter));
    mxSetField(mxPtr, 0, "eps_abs",                   mxCreateDoubleScalar(settings->eps_abs));
    mxSetField(mxPtr, 0, "eps_abs_in",                mxCreateDoubleScalar(settings->eps_abs_in));
    mxSetField(mxPtr, 0, "eps_prim_inf",              mxCreateDoubleScalar(settings->eps_prim_inf));
    mxSetField(mxPtr, 0, "eps_dual_inf",              mxCreateDoubleScalar(settings->eps_dual_inf));
    mxSetField(mxPtr, 0, "rho",                       mxCreateDoubleScalar(settings->rho));
    mxSetField(mxPtr, 0, "theta",                     mxCreateDoubleScalar(settings->theta));
    mxSetField(mxPtr, 0, "delta",                     mxCreateDoubleScalar(settings->delta));
    mxSetField(mxPtr, 0, "mu_min",                    mxCreateDoubleScalar(settings->mu_min));
    mxSetField(mxPtr, 0, "proximal",                  mxCreateDoubleScalar(settings->proximal));
    mxSetField(mxPtr, 0, "sigma_init",                mxCreateDoubleScalar(settings->sigma_init));
    mxSetField(mxPtr, 0, "sigma_upd",                 mxCreateDoubleScalar(settings->sigma_upd));
    mxSetField(mxPtr, 0, "sigma_min",                 mxCreateDoubleScalar(settings->sigma_min));
    mxSetField(mxPtr, 0, "scaling",                   mxCreateDoubleScalar(settings->scaling));
    mxSetField(mxPtr, 0, "verbose",                   mxCreateDoubleScalar(settings->verbose));
    mxSetField(mxPtr, 0, "print_interval",            mxCreateDoubleScalar(settings->print_interval));
    mxSetField(mxPtr, 0, "reset_newton_iter",         mxCreateDoubleScalar(settings->reset_newton_iter));
    mxSetField(mxPtr, 0, "max_time",                  mxCreateDoubleScalar(settings->max_time));

    return mxPtr;
}


void copyMxStructToSettings(const mxArray* mxPtr, QPDOSettings* settings){

    // this function assumes that only a complete and validated structure
    // will be passed.  matlab mex-side function is responsible for checking
    // structure validity

    // map the QPDO_SETTINGS fields one at a time into mxArrays
    // matlab handles everything as a double
    settings->max_iter                  = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "max_iter"));
    settings->inner_max_iter            = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "inner_max_iter"));
    settings->eps_abs                   = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs"));
    settings->eps_abs_in                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs_in"));
    settings->eps_prim_inf              = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_prim_inf"));
    settings->eps_dual_inf              = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
    settings->rho                       = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "rho"));
    settings->theta                     = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "theta"));
    settings->delta                     = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
    settings->mu_min                    = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "mu_min"));
    settings->proximal                  = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "proximal"));
    settings->sigma_init                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "sigma_init"));
    settings->sigma_upd                 = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "sigma_upd"));
    settings->sigma_min                 = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "sigma_min"));
    settings->scaling                   = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "scaling"));
    settings->verbose                   = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
    settings->print_interval            = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "print_interval"));
    settings->reset_newton_iter         = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "reset_newton_iter"));
    settings->max_time                  = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "max_time"));

}
