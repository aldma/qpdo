#include "iteration.h"
#include "lin_alg.h"
#include "cholmod_interface.h"
#include "newton.h"
#include "linesearch.h"

/*************************************************************************
 * UPDATE_ITERATE
 * search direction, exact linesearch, update
 *************************************************************************/
void update_iterate(QPDOWorkspace *work, cholmod_common *c) {
    size_t n = work->data->n;
    size_t m = work->data->m;
    // search direction (semismooth Newton's)
    newton_direction(work, c);
    // line search (exact, sorting)
    exact_linesearch(work, c);
    // update x and y
    vec_add_scaled(work->x, work->dx, work->x, work->tau, n);
    vec_add_scaled(work->y, work->dy, work->y, work->tau, m);
    // update Qx, Ax, Aty
    vec_add_scaled(work->Qx, work->Qdx, work->Qx, work->tau, n);
    vec_add_scaled(work->Ax, work->Adx, work->Ax, work->tau, m);
    vec_add_scaled(work->Aty, work->Atdy, work->Aty, work->tau, n);
}

/*************************************************************************
 * COMPUTE_OUTER_RESIDUALS
 *************************************************************************/
void compute_outer_residuals(QPDOWorkspace *work, cholmod_common *c) {
    // dimensions
    size_t n = work->data->n;
    size_t m = work->data->m;
    size_t i;
    // primal residuals : feasibility
    // res_prim <- Ax - max( l, min( Ax + y, u ))
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->E, work->y, work->temp_m, m);
        vec_self_mult_scalar(work->temp_m, work->scaling->cinv, m);
        vec_ew_prod(work->scaling->E, work->temp_m, work->temp_m, m);
        vec_add_scaled(work->Ax, work->temp_m, work->temp_m, 1, m); // Ax + y
    } else {
        vec_add_scaled(work->Ax, work->y, work->temp_m, 1, m); // Ax + y
    }
    vec_ew_mid_vec(work->temp_m, work->data->l, work->data->u, work->temp_m, m);
    vec_add_scaled(work->Ax, work->temp_m, work->res_prim, -1, m);
    // dual residuals : optimality
    // res_dual = Q*x + q + A'*y
    // df <- { Q * x + q + sigma * x , if proximal
    //       { Q * x + q             , otherwise
    // NB: temporary `df`
    vec_add_scaled(work->Qx, work->data->q, work->df, 1, n);
    // res_dual
    if (work->settings->proximal) {
        vec_add_scaled(work->df, work->x, work->res_dual, - work->sigma, n);
        vec_add_scaled(work->res_dual, work->Aty, work->res_dual, 1, n);
    } else {
        vec_add_scaled(work->df, work->Aty, work->res_dual, 1, n);
    }
}

/*************************************************************************
 * COMPUTE_INNER_RESIDUALS
 *************************************************************************/
void compute_inner_residuals(QPDOWorkspace *work, cholmod_common *c) {
    // dimensions
    size_t n = work->data->n;
    size_t m = work->data->m;
    size_t i;
    // primal residuals : feasibility
    // res_prim_in = A * x + mu .* (ybar - y) - proj(A * x + mu .* (ybar - y/2))
    // w <- Ax + mu .* (ybar - y/2)
    for (i = 0; i < m; i++) {
        work->w[i] = work->Ax[i] + work->mu[i] * ( work->ybar[i] - 0.5*work->y[i] );
    }
    // zin = max( l, min( w, u ))
    vec_ew_mid_vec(work->w, work->data->l, work->data->u, work->temp_m, m);
    // res_prim_in <- Ax + mu .* (ybar - y) - zin
    for (i = 0; i < m; i++) {
        work->res_prim_in[i] = work->Ax[i] + work->mu[i] * ( work->ybar[i] - work->y[i] ) - work->temp_m[i];
    }
    // dual residuals : optimality
    // res_dual_in = { Q * x + q + A' * y + sigma * (x - xbar) , if proximal
    //               { Q * x + q + A' * y                      , otherwise
    // NB: here we assume `df` already contains `Qx + q`
    // df <- { Q*x + q + sigma * (x-xbar) , if proximal
    //       { Q*x + q                    , otherwise
    if (work->settings->proximal) {
        vec_add_scaled(work->df, work->xbar, work->df, - work->sigma, n);
    }
    // res_dual_in <- df + Aty
    vec_add_scaled(work->df, work->Aty, work->res_dual_in, 1, n);
}

/*************************************************************************
 * INITIALIZE_MU
 *************************************************************************/
void initialize_mu(QPDOWorkspace *work, cholmod_common *c) {
    size_t n = work->data->n;
    size_t m = work->data->m;
    // Compute initial dual regularization/penalty parameter ´mu´
    c_float f = 0.5*vec_prod(work->x, work->Qx, n) + vec_prod(work->data->q, work->x, n); // 1/2 x'*Q*x + x'*q
    vec_ew_mid_vec(work->Ax, work->data->l, work->data->u, work->temp_m, m); // z : max(l, min(A*x, u))
    vec_add_scaled(work->Ax, work->temp_m, work->temp_m, -1, m); // A*x - z
    // balancing penalty parameter from [Birgin & Martinez, 2014, §12.4]
    size_t i;
    for (i = 0; i < m; i++) {
        work->mu[i] = c_max( 1e-3, c_min( 1e3, 0.1 * c_max( 1, 0.5 * work->temp_m[i] * work->temp_m[i] ) / c_max( 1, c_absval(f) ) ) );
    }
    // Set fields related to mu
    // sqrt_mu = 1 ./ sqrt(mu)
    vec_ew_sqrt(work->mu, work->sqrt_mu, m);
    vec_ew_recipr(work->sqrt_mu, work->sqrt_mu, m);
    // At_scale
    c_float *At_scalex = work->chol->At_scale->x;
    prea_vec_copy(work->sqrt_mu, At_scalex, m);
    if (work->chol->At_sqrt_mu) CHOLMOD(free_sparse)(&work->chol->At_sqrt_mu, c);
    work->chol->At_sqrt_mu = CHOLMOD(transpose)(work->data->A, 1, c);
    CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_mu, c);
    // sqrt_mu_min = 1 / sqrt( mu_min )
    work->sqrt_mu_min = 1 / c_sqrt( work->settings->mu_min );
}

/*************************************************************************
 * UPDATE_MU
 *************************************************************************/
void update_mu(QPDOWorkspace* work, cholmod_common *c) {
    work->n_mu_changed = 0;
    c_float *At_scalex = work->chol->At_scale->x;
    c_float res_prim_norm = vec_norm_inf(work->res_prim, work->data->m);
    c_float mu_new, mu_factor;
    c_int *mu_changed = work->chol->enter;
    size_t k;
    for (k = 0; k < work->data->m; k++) {
        if (c_absval(work->res_prim[k]) > c_max(work->settings->eps_abs, work->settings->theta * c_absval(work->res_prim_old[k]))) {
            mu_factor = 1.0 / c_min( 1.0, work->settings->delta * res_prim_norm / c_absval(work->res_prim[k]) );
            mu_new = work->mu[k] / mu_factor;
            if (mu_new >= work->settings->mu_min) {
                if (work->mu[k] != mu_new) {
                    mu_changed[work->n_mu_changed] = (c_int)k;
                    work->n_mu_changed++;
                }
                work->mu[k] = mu_new;
                mu_factor = c_sqrt( mu_factor );
                work->sqrt_mu[k] = mu_factor * work->sqrt_mu[k];
                At_scalex[k] = mu_factor;
            } else {
                if (work->mu[k] != work->settings->mu_min) {
                    mu_changed[work->n_mu_changed] = (c_int)k;
                    work->n_mu_changed++;
                }
                work->mu[k] = work->settings->mu_min;
                At_scalex[k] = work->sqrt_mu_min / work->sqrt_mu[k];
                work->sqrt_mu[k] = work->sqrt_mu_min;
            }
        } else {
            At_scalex[k] = 1.0; // no update
        }
    }
    CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_mu, c);
    if ((work->settings->proximal && work->sigma > work->settings->sigma_min) || (work->n_mu_changed > 0.25*MAX_RANK_UPDATE)) {
        work->chol->reset_newton = TRUE;
    } else if (work->n_mu_changed == 0){
        // do nothing
    } else {
        ldlupdate_mu_changed(work, c);
    }
}

/*************************************************************************
 * UPDATE_SIGMA
 *************************************************************************/
void update_sigma(QPDOWorkspace *work) {
    if (work->sigma > work->settings->sigma_min) {
        c_float sigma_old = work->sigma;
        work->sigma = c_max( work->sigma * work->settings->sigma_upd, work->settings->sigma_min );
        work->chol->reset_newton = TRUE;
        vec_add_scaled(work->Qx, work->x, work->Qx, work->sigma - sigma_old, work->data->n);
    }
}

/*************************************************************************
 * COMPUTE_OBJECTIVE
 *************************************************************************/
c_float compute_objective(QPDOWorkspace *work) {
    c_float objective = 0;
    size_t n = work->data->n;
    size_t i = 0;
    if (work->settings->proximal) {
        if(n >= 4) {
            for (; i <= n-4; i+=4) {
                objective +=  (0.5*(work->Qx[i]   - work->x[i]   * work->sigma) + work->data->q[i]  )*work->x[i]
                            + (0.5*(work->Qx[i+1] - work->x[i+1] * work->sigma) + work->data->q[i+1])*work->x[i+1]
                            + (0.5*(work->Qx[i+2] - work->x[i+2] * work->sigma) + work->data->q[i+2])*work->x[i+2]
                            + (0.5*(work->Qx[i+3] - work->x[i+3] * work->sigma) + work->data->q[i+3])*work->x[i+3];
            } /* cache-friendly implementation */
        }
        for (; i < n; i++) {
            objective += (0.5*(work->Qx[i] - work->sigma * work->x[i])+ work->data->q[i])*work->x[i];
        }
    } else {
        if(n >= 4) {
            for (; i <= n-4; i+=4) {
                objective +=  (0.5*work->Qx[i]   + work->data->q[i]  )*work->x[i]
                            + (0.5*work->Qx[i+1] + work->data->q[i+1])*work->x[i+1]
                            + (0.5*work->Qx[i+2] + work->data->q[i+2])*work->x[i+2]
                            + (0.5*work->Qx[i+3] + work->data->q[i+3])*work->x[i+3];
            } /* cache-friendly implementation */
        }
        for (; i < n; i++) {
            objective += (0.5*work->Qx[i] + work->data->q[i])*work->x[i];
        }/* include remaining terms at the end */
    }
    // scaling
    if (work->settings->scaling) {
        objective *= work->scaling->cinv;
    }
    // constant term
    objective += work->data->c;
    return objective;
}
