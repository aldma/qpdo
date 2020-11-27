#include "termination.h"
#include "lin_alg.h"
#include "constants.h"
#include "global_opts.h"
#include "util.h"
#include "iteration.h"

/*************************************************************************
 * CHECK_OUTER_OPTIMALITY
 *************************************************************************/
c_int check_outer_optimality(QPDOWorkspace *work) {
    // check divergence
    if ((work->info->res_prim_norm > QPDO_INFTY) || (work->info->res_dual_norm > QPDO_INFTY)) {
        update_status(work->info, QPDO_NON_CVX);
        return 1;
    }
    // check optimality
    if ((work->info->res_prim_norm <= work->settings->eps_abs) && (work->info->res_dual_norm <= work->settings->eps_abs)) {
        update_status(work->info, QPDO_SOLVED);
        return 1;
    }
    return 0;
}

/*************************************************************************
 * CHECK_INNER_OPTIMALITY
 *************************************************************************/
c_int check_inner_optimality(QPDOWorkspace *work) {
    return (work->info->res_prim_in_norm <= work->eps_in) && (work->info->res_dual_in_norm <= work->eps_in);
}

/*************************************************************************
 * COMPUTE_OUTER_RESIDUALS_NORM
 *************************************************************************/
void compute_outer_residuals_norm(QPDOWorkspace *work) {
    size_t n = work->data->n;
    size_t m = work->data->m;
    if (work->settings->scaling) {
        // primal
        vec_ew_prod(work->scaling->Einv, work->res_prim, work->temp_m, m);
        work->info->res_prim_norm = vec_norm_inf(work->temp_m, m);
        // dual
        vec_ew_prod(work->scaling->Dinv, work->res_dual, work->temp_n, n);
        work->info->res_dual_norm = vec_norm_inf(work->temp_n, n);
        work->info->res_dual_norm *= work->scaling->cinv;
    } else {
        // primal
        work->info->res_prim_norm = vec_norm_inf(work->res_prim, m);
        // dual
        work->info->res_dual_norm = vec_norm_inf(work->res_dual, n);
    }
    return;
}

/*************************************************************************
 * COMPUTE_INNER_RESIDUALS_NORM
 *************************************************************************/
void compute_inner_residuals_norm(QPDOWorkspace *work, cholmod_common *c) {
    size_t n = work->data->n;
    size_t m = work->data->m;
    // primal
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Einv, work->res_prim_in, work->temp_m, m);
        work->info->res_prim_in_norm = vec_norm_inf(work->temp_m, m);
    } else {
        work->info->res_prim_in_norm = vec_norm_inf(work->res_prim_in, m);
    }
    // dual
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Dinv, work->res_dual_in, work->temp_n, n);
        work->info->res_dual_in_norm = vec_norm_inf(work->temp_n, n);
        work->info->res_dual_in_norm *= work->scaling->cinv;
    } else {
        work->info->res_dual_in_norm = vec_norm_inf(work->res_dual_in, n);
    }
    return;
}

/*************************************************************************
 * STORE_SOLUTION
 *************************************************************************/
void store_solution(QPDOWorkspace *work) {
    if (work->settings->scaling) {
        vec_ew_prod(work->x, work->scaling->D, work->solution->x, work->data->n);
        vec_self_mult_scalar(work->y, work->scaling->cinv, work->data->m);
        vec_ew_prod(work->y, work->scaling->E, work->solution->y, work->data->m);
    } else {
        prea_vec_copy(work->x, work->solution->x, work->data->n);
        prea_vec_copy(work->y, work->solution->y, work->data->m);
    }
    work->info->objective = compute_objective( work );
}
