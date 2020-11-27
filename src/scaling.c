# ifdef __cplusplus
extern "C" {
# endif

#include "scaling.h"
#include "cholmod.h"
#include <stdio.h>

/*************************************************************************
 * LIMIT_SCALING
 * set values lower than threshold to 1
 *************************************************************************/
void limit_scaling(c_float *D, size_t n) {
  size_t i;
  for (i = 0; i < n; i++) {
    D[i] = D[i] < MIN_SCALING ? 1.0 : D[i];
  }
}

/*************************************************************************
 * CHECK_OUTER_OPTIMALITY
 * Ruiz equilibration of constraint matrix + cost scaling
 *************************************************************************/
void scale_data(QPDOWorkspace *work) {
    cholmod_common c;
    CHOLMOD(start)(&c);
    cholmod_set_settings(&c);

    size_t n = work->data->n;
    size_t m = work->data->m;
    vec_set_scalar(work->scaling->D, 1, n);
    vec_set_scalar(work->scaling->E, 1, m);

    c_int i;
    // Ruiz on constraint matrix A
    for (i = 0; i < work->settings->scaling; i++) {

        // Set D_temp = vecnorm(A,inf,1) (cols) and E_temp = vecnorm(A,inf,2) (rows)
        mat_inf_norm_cols(work->data->A, work->D_temp);
        mat_inf_norm_rows(work->data->A, work->E_temp);

        // Set to 1 values with 0 norms (avoid crazy scaling)
        limit_scaling(work->D_temp, n);
        limit_scaling(work->E_temp, m);

        // Take square root of norms
        vec_ew_sqrt(work->D_temp, work->D_temp, n);
        vec_ew_sqrt(work->E_temp, work->E_temp, m);

        // 1./D and 1./E
        vec_ew_recipr(work->D_temp, work->D_temp, n);
        vec_ew_recipr(work->E_temp, work->E_temp, m);

        // Equilibrate matrix A
        // A <- E A D
        CHOLMOD(scale)(work->chol->E_temp, CHOLMOD_ROW, work->data->A, &c);
        CHOLMOD(scale)(work->chol->D_temp, CHOLMOD_COL, work->data->A, &c);

        // Update equilibration matrices D and E
        vec_ew_prod(work->scaling->D, work->D_temp, work->scaling->D, n);
        vec_ew_prod(work->scaling->E, work->E_temp, work->scaling->E, m);

    }

    // Equilibrate matrix Q and vector q
    // Q <- D Q D, q <- D q
    prea_vec_copy(work->scaling->D, work->D_temp, n);
    CHOLMOD(scale)(work->chol->D_temp, CHOLMOD_SYM, work->data->Q, &c);
    vec_ew_prod(work->scaling->D, work->data->q, work->data->q, n);

    // Cost scaling
    vec_add_scaled(work->Qx, work->data->q, work->temp_n, 1, n);
    work->scaling->c = 1/c_max(1.0, vec_norm_inf(work->temp_n, n));
    vec_self_mult_scalar(work->data->q, work->scaling->c, n);
    cholmod_dense *scalar = CHOLMOD(ones)(1,1,CHOLMOD_REAL, &c);
    c_float *scalarx = scalar->x;
    scalarx[0] = work->scaling->c;
    CHOLMOD(scale)(scalar, CHOLMOD_SCALAR, work->data->Q, &c);
    CHOLMOD(free_dense)(&scalar, &c);

    // Store cinv, Dinv, Einv
    vec_ew_recipr(work->scaling->D, work->scaling->Dinv, n);
    vec_ew_recipr(work->scaling->E, work->scaling->Einv, m);
    work->scaling->cinv = (c_float) 1.0/work->scaling->c;

    // Scale problem vectors l, u
    vec_ew_prod(work->scaling->E, work->data->l, work->data->l, m);
    vec_ew_prod(work->scaling->E, work->data->u, work->data->u, m);

    CHOLMOD(finish)(&c);
}


# ifdef __cplusplus
}
# endif
