#ifndef QPDO_TYPES_H
# define QPDO_TYPES_H

# ifdef __cplusplus
extern "C" {
# endif

#include "global_opts.h"
#include "cholmod.h"

/*************************************************************************
 * LINESEARCH ARRAY
 *************************************************************************/
typedef struct array_element  {
  c_float x; // value of the element
  size_t  i; // index
} array_element;


/*************************************************************************
 * QPDO TYPES
 *************************************************************************/

/*************************************************************************
 * SOLUTION
 *************************************************************************/
typedef struct {
  c_float *x; // primal solution
  c_float *y; // dual solution
} QPDOSolution;

/*************************************************************************
 * TIMER
 *************************************************************************/
typedef struct QPDO_TIMER QPDOTimer;

/*************************************************************************
 * SCALING
 *************************************************************************/
typedef struct {
  c_float *D;    // primal variable scaling
  c_float *Dinv; // primal variable rescaling
  c_float *E;    // dual variable scaling
  c_float *Einv; // dual variable rescaling
  c_float  c;    // objective scaling
  c_float  cinv; // objective rescaling
} QPDOScaling;


/*************************************************************************
 * INFORMATION
 *************************************************************************/
typedef struct {
  c_int   iterations;     // number of iterations taken
  c_int   oterations;     // number of outer iterations (i.e. dual updates)
  char    status[32];     // status string, e.g. 'solved'
  c_int   status_val;     // status as c_int, defined in constants.h

  c_float res_prim_norm;    // norm of primal residual
  c_float res_dual_norm;    // norm of dual residual
  c_float res_prim_in_norm; // norm of inner primal residual
  c_float res_dual_in_norm; // norm of inner dual residual

  c_float objective;      // objective function value

  #ifdef PROFILING
  c_float setup_time;    // time taken for setup phase (seconds)
  c_float solve_time;    // time taken for solve phase (seconds)
  c_float run_time;      // total time (seconds)
  #endif

} QPDOInfo;

/**********************************
* Main structures and Data Types *
**********************************/

/**
 * Data structure
 */
typedef struct {
  size_t          n;    // number of variables n
  size_t          m;    // number of constraints m
  cholmod_sparse *Q;    // sparse quadratic part of the cost Q (size n x n)
  cholmod_sparse *A;    // sparse linear constraints matrix A (size m x n)
  c_float        *q;    // dense array for linear part of cost function (size n)
  c_float         c;    // constant part of cost
  c_float        *l;    // dense array for lower bounds (size m)
  c_float        *u;    // dense array for upper bounds (size m)
} QPDOData;


/**
 * Settings struct
 */
typedef struct {
  c_float max_time;                 // time limit
  c_int   max_iter;                 // maximum number of iterations
  c_int   inner_max_iter;           // maximum number of inner iterations
  c_float eps_abs;                  // absolute tolerance
  c_float eps_abs_in;               // initial inner absolute tolerance
  c_float eps_prim_inf;             // primal infeasibility tolerance
  c_float eps_dual_inf;             // dual infeasibility tolerance
  c_float rho;                      // inner tolerance shrink factor
  c_float theta;                    // penalty update check factor
  c_float delta;                    // penalty update shrink factor
  c_float mu_min;                   // minimum dual/penalty parameter
  c_int   proximal;                 // boolean, use primal regularization?
  c_float sigma_init;               // initial proximal penalty parameter
  c_float sigma_upd;                // proximal penalty update factor
  c_float sigma_min;                // proximal penalty parameter cap
  c_int   scaling;                  // scaling iterations, if 0 then scaling is disabled
  c_int   verbose;                  // boolean, write out progress
  c_int   print_interval;           // frequency of printing
  c_int   reset_newton_iter;        // frequency of performing a complete Cholesky factorization
} QPDOSettings;

/*************************************************************************
 * LINEAR SOLVER (Cholmod)
 *************************************************************************/
typedef struct {
  cholmod_factor *LD;             // LD factor (part of LDL' factorization)
  cholmod_factor *LD_Q;           // LD factor of Q (useful in computing dual objective)
  cholmod_dense *E_temp;          // temporary constraints scaling vectors
  cholmod_dense *D_temp;          // temporary primal variable scaling vectors
  cholmod_dense *linsys_rhs;      // rhs
  cholmod_dense *dx;              // primal search direction
  cholmod_dense *dy;              // dual search direction
  cholmod_dense *Adx;             // A * d
  cholmod_dense *Qdx;             // Q * d
  cholmod_dense *Atdy;            // A' * dy
  c_int reset_newton;             // boolean, after mu is updated perform a new factorization
  c_int *active_set;              // index set of active constraints
  c_int *active_set_old;          // index set of active constraints in the previous iteration
  c_int n_active_set;             // number of active constraints
  c_int *enter;                   // index set of entering constraints
  c_int n_enter;                  // number of entering constraints
  c_int *leave;                   // index set of leaving constraints
  c_int n_leave;                  // number of leaving constraints
  cholmod_dense *At_scale;        // running vector of sqrt(mu), used to scale At_sqrt_mu
  cholmod_sparse *At_sqrt_mu;     // A' / sqrt(mu)
} QPDOCholmod;

/*************************************************************************
 * WORKSPACE
 *************************************************************************/
typedef struct {
  QPDOData *data; // problem data to work on (possibly scaled)

  // iterate
  c_float *x;        // primal variable
  c_float *y;        // dual variable
  c_float *Ax;       // A * x
  c_float *Qx;       // Q * x
  c_float *Aty;      // A' * y
  c_int initialized; // flag whether the iterates were initialized or not

  // temporaries
  c_float *temp_m;          // m-vector
  c_float *temp_n;          // n-vector
  c_float *temp_2m;         // 2m-vector

  // parameters
  c_float *mu;              // penalty vector
  c_float *sqrt_mu;         // elementwise sqrt(mu)
  c_float sqrt_mu_min;      // sqrt(mu_min)
  c_float sqrt_delta;       // sqrt(penalty update factor)
  c_int   n_mu_changed;     // number of mu-components that changed in an outer iteration
  c_float sigma;            // proximal penalty factor
  c_int   sigma_mined;      // flag to indicate whether sigma has been minimized when the primal residual was low
  c_float norm_q;

  // base point and step
  c_float *xbar;          // primal estimate
  c_float *ybar;          // dual estimate
  c_float *dx;            // primal search direction / step
  c_float *dy;            // dual search direction / step
  c_float tau;          // stepsize
  c_float *Qdx;         // Q * dx
  c_float *Adx;         // A * dx
  c_float *Atdy;        // A' * dy

  c_float *w;              // Ax + mu .* y
  c_float *z;              // projection of w onto the constraint set [l, u]
  c_float *df;             // gradient of the primal objective (+ proximal term)
  c_float *res_prim;       // primal residual
  c_float *res_dual;       // dual residual
  c_float *res_prim_old;   // old primal residual
  c_float *res_prim_in;    // inner primal residual
  c_float *res_dual_in;    // inner dual residual
  c_float *linsys_rhs;     // linear system rhs

  c_float res_prim_norm_old;
  c_float res_dual_norm_old;

  c_float ls_eta;         // linesearch parameter
  c_float ls_beta;        // linesearch parameter
  c_float *ls_delta;      // linesearch parameter
  c_float *ls_alpha;      // linesearch parameter
  array_element *ls_taus; // alpha ./ delta
  c_int   *ls_idx_L;      // index set L (where s>0)
  c_int   *ls_idx_P;      // index set P (where delta>0)
  c_int   *ls_idx_J;      // index set J (L xor P)

  c_float eps_prim;      // primal tolerance
  c_float eps_dual;      // dual tolerance
  c_float eps_prim_in;   // inner primal tolerance
  c_float eps_dual_in;   // inner dual tolerance
  c_float eps_in;        // inner absolute tolerance

  c_float *D_temp;      // temporary primal variable scaling vectors
  c_float *E_temp;      // temporary constraints scaling vectors

  QPDOCholmod  *chol;     // cholmod variables
  QPDOSettings *settings; // problem settings
  QPDOScaling  *scaling;  // scaling vectors
  QPDOSolution *solution; // problem solution
  QPDOInfo     *info;     // solver information

  # ifdef PROFILING
    QPDOTimer *timer;     // timer object
  # endif

} QPDOWorkspace;


# ifdef __cplusplus
}
# endif

#endif
