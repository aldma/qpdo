#ifndef CHOLMOD_INTERFACE_H
#define CHOLMOD_INTERFACE_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "cholmod.h"
#include "global_opts.h"
#include "constants.h"
#include "types.h"


/**
 * Matrix-vector multiplication.
 */
void mat_vec(cholmod_sparse *A,
             cholmod_dense  *x,
             cholmod_dense  *y,
             cholmod_common *c);
/**
 * Matrix-transpose-vector multiplication.
 */
void mat_tpose_vec(cholmod_sparse *A, cholmod_dense  *x, cholmod_dense  *y, cholmod_common *c);
/**
 * Infinity norm of each matrix column, E_i = \| M{(:,i)} \|_\infty.
 */
void mat_inf_norm_cols(cholmod_sparse *M, c_float   *E);

/**
 * Infinity norm of each matrix row, E_i = \| M{(i,:)} \|_\infty.
 */
void mat_inf_norm_rows(cholmod_sparse *M, c_float   *E);

/**
 * Calculate LDL^T factorization of a matrix M.
 */
void ldlchol(cholmod_sparse *M, QPDOWorkspace *work, cholmod_common *c);

/**
 * Calculate LDL^T factorization of Q + A{(act,:)}^T * \Mu{(act,act)}^{-1} * A{(act,:)}, with \Mu=diag(\mu) and act the set of active constraints.
 */
void ldlcholQAtmuA(QPDOWorkspace *work, cholmod_common *c);

/**
 * Update the LDL^T factorization given a set of entering constraints.
 */
void ldlupdate_enter_constraints(QPDOWorkspace *work, cholmod_common *c);

/**
 * Downdate the LDL^T factorization given a set of leaving constraints.
 */
void ldldowndate_leave_constraints(QPDOWorkspace *work, cholmod_common *c);

/**
 * Update the LDL^T factorization given a set of indexes where mu has been updated.
 */
void ldlupdate_mu_changed(QPDOWorkspace *work, cholmod_common *c);


/**
 * Solve the linear system LDL^T d = rhs.
 */
void ldlsolveLD_rhs(QPDOWorkspace *work, cholmod_common *c);

/**
 * Cholmod settings and memory allocation.
 */
void cholmod_set_settings(cholmod_common *c);

# ifdef __cplusplus
}
# endif

#endif
