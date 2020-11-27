#include "cholmod_interface.h"
#include "lin_alg.h"
#include <stdio.h>

/*************************************************************************
 * LDLCHOL                                                               *
 *************************************************************************/
void ldlchol(cholmod_sparse *M, QPDOWorkspace *work, cholmod_common *c) {
    if (work->chol->LD) CHOLMOD(free_factor)(&work->chol->LD, c);
    work->chol->LD = CHOLMOD(analyze) (M, c) ;
    if (work->settings->proximal) {
        double beta [2] = {work->sigma,0};
        CHOLMOD(factorize_p) (M, beta, NULL, 0, work->chol->LD, c);
    } else {
        CHOLMOD(factorize) (M, work->chol->LD, c);
    }
    // If integers are used, supernodal might fail, so check for this and
    // switch to simplicial if necessary.
    #ifndef DLONG
        if ((c)->status != CHOLMOD_OK) {
            (c)->supernodal = CHOLMOD_SIMPLICIAL;
            if (work->settings->proximal) {
                double beta [2] = {work->sigma,0};
                CHOLMOD(factorize_p) (M, beta, NULL, 0, work->chol->LD, c);
            } else {
                CHOLMOD(factorize) (M, work->chol->LD, c);
            }
        }
    #endif
}

/*************************************************************************
 * LDLCHOLQATMUA                                                         *
 *************************************************************************/
void ldlcholQAtmuA(QPDOWorkspace *work, cholmod_common *c) {
    cholmod_sparse *AtmuA;
    cholmod_sparse *QAtmuA;
    size_t n_active = 0;
    for (size_t i = 0; i < work->data->m; i++) {
        if (work->chol->active_set[i]){
            work->chol->enter[n_active] = (c_int) i;
            n_active++;
        }
    }
    AtmuA = CHOLMOD(aat)(work->chol->At_sqrt_mu, work->chol->enter, n_active, TRUE, c);
    double one [2] = {1,0};
    QAtmuA = CHOLMOD(add)(work->data->Q, AtmuA, one, one, TRUE, FALSE, c);
    QAtmuA->stype = work->data->Q->stype;
    ldlchol(QAtmuA, work, c);
    CHOLMOD(free_sparse)(&AtmuA, c);
    CHOLMOD(free_sparse)(&QAtmuA, c);
}

/*************************************************************************
 * LDLUPDATE_ENTERING_CONSTRAINTS                                        *
 *************************************************************************/
void ldlupdate_enter_constraints(QPDOWorkspace *work, cholmod_common *c) {
    cholmod_sparse *Ae;
    Ae = CHOLMOD(submatrix)(work->chol->At_sqrt_mu, NULL, -1, work->chol->enter, work->chol->n_enter, TRUE, TRUE, c);
    CHOLMOD(updown)(TRUE, Ae, work->chol->LD, c);
    CHOLMOD(free_sparse)(&Ae, c);
}

/*************************************************************************
 * LDLDOWNDATE_LEAVING_CONSTRAINTS                                       *
 *************************************************************************/
void ldldowndate_leave_constraints(QPDOWorkspace *work, cholmod_common *c) {
    cholmod_sparse *Al;
    Al = CHOLMOD(submatrix)(work->chol->At_sqrt_mu, NULL, -1, work->chol->leave, work->chol->n_leave, TRUE, TRUE, c);
    CHOLMOD(updown)(FALSE, Al, work->chol->LD, c);
    CHOLMOD(free_sparse)(&Al, c);
}

/*************************************************************************
 * LDLUPDATE_MU_CHANGED                                                  *
 *************************************************************************/
void ldlupdate_mu_changed(QPDOWorkspace *work, cholmod_common *c) {
    cholmod_sparse *Ae;
    c_float *At_scalex = work->chol->At_scale->x;
    c_int *mu_changed = work->chol->enter;
    size_t k, n_mu_changed = (size_t) work->n_mu_changed;
    for (k=0; k < n_mu_changed; k++) {
        At_scalex[mu_changed[k]] = c_sqrt( 1 - 1/(At_scalex[mu_changed[k]] * At_scalex[mu_changed[k]]) );
    }
    CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_mu, c);
    Ae = CHOLMOD(submatrix)(work->chol->At_sqrt_mu, NULL, -1, mu_changed, work->n_mu_changed, TRUE, TRUE, c);
    for (k=0; k < work->data->m; k++) {
        At_scalex[k] = 1.0 / At_scalex[k];
    }
    CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_mu, c);
    CHOLMOD(updown)(TRUE, Ae, work->chol->LD, c);
    CHOLMOD(free_sparse)(&Ae, c);
}

/*************************************************************************
 * LDLSOLVELD_RHS                                                        *
 *************************************************************************/
void ldlsolveLD_rhs(QPDOWorkspace *work, cholmod_common *c) {
    if (work->chol->dx) CHOLMOD(free_dense)(&work->chol->dx, c);
    work->chol->dx = CHOLMOD(solve) (CHOLMOD_LDLt, work->chol->LD, work->chol->linsys_rhs, c);
    work->dx = work->chol->dx->x;
}

/*************************************************************************
 * CHOLMOD_SET_SETTINGS                                                  *
 *************************************************************************/
void cholmod_set_settings(cholmod_common *c) {
    SuiteSparse_config.malloc_func = c_malloc;
    SuiteSparse_config.calloc_func = c_calloc;
    SuiteSparse_config.free_func = c_free;
    SuiteSparse_config.realloc_func = c_realloc;
    c->final_asis = FALSE ;
    c->final_super = FALSE ;
    c->final_ll = FALSE ;
    c->final_pack = TRUE ;
    c->final_monotonic = TRUE ;
    c->final_resymbol = TRUE ;
    c->quick_return_if_not_posdef = TRUE;
    c->nmethods = 1 ;
    c->method [0].ordering = CHOLMOD_NATURAL ;
    c->postorder = FALSE ;
    c->useGPU = FALSE ;
}

/*************************************************************************
 * MATRIX OPERATIONS                                                     *
 *************************************************************************/

/*************************************************************************
 * MAT_VEC                                                               *
 *************************************************************************/
void mat_vec(cholmod_sparse *A, cholmod_dense *x, cholmod_dense *y, cholmod_common *c) {
    double one [2] = {1,0};
    double zero [2] = {0,0};
    if (x!=y) {
        CHOLMOD(sdmult)(A, 0, one, zero, x, y, c);
    } else {
        cholmod_dense* x2 = CHOLMOD(copy_dense)(x, c);
        CHOLMOD(sdmult)(A, 0, one, zero, x2, y, c);
        CHOLMOD(free_dense)(&x2, c);
    }
}

/*************************************************************************
 * MAT_TPOSE_VEC                                                         *
 *************************************************************************/
void mat_tpose_vec(cholmod_sparse *A, cholmod_dense *x, cholmod_dense *y, cholmod_common *c) {
    double one [2] = {1,0};
    double zero [2] = {0,0};
    if (x!=y) {
        CHOLMOD(sdmult)(A, 1, one, zero, x, y, c);
    } else {
        cholmod_dense* x2 = CHOLMOD(copy_dense)(x, c);
        CHOLMOD(sdmult)(A, 1, one, zero, x2, y, c);
        CHOLMOD(free_dense)(&x2, c);
    }
}

/*************************************************************************
 * MAT_INF_NORM_COLS                                                     *
 *************************************************************************/
void mat_inf_norm_cols(cholmod_sparse *M, c_float *E) {
    size_t j;
    c_int k;
    c_float *Mx = M->x;
    c_int *Mp = M->p;
    // Initialize zero max elements
    for (j = 0; j < M->ncol; j++) {
        E[j] = 0.;
    }
    // Compute maximum across columns
    for (j = 0; j < M->ncol; j++) {
        for (k = Mp[j]; k < Mp[j + 1]; k++) {
            E[j] = c_max(c_absval(Mx[k]), E[j]);
        }
    }
}

/*************************************************************************
 * MAT_INF_NORM_ROWS                                                     *
 *************************************************************************/
void mat_inf_norm_rows(cholmod_sparse *M, c_float *E) {
    size_t j;
    c_int i, k;
    c_int *Mp = M->p;
    c_int *Mi = M->i;
    c_float *Mx = M->x;
    // Initialize zero max elements
    for (j = 0; j < M->nrow; j++) {
        E[j] = 0.;
    }
    // Compute maximum across rows
    for (j = 0; j < M->ncol; j++) {
        for (k = Mp[j]; k < Mp[j + 1]; k++) {
            i    = Mi[k];
            E[i] = c_max(c_absval(Mx[k]), E[i]);
        }
    }
}
