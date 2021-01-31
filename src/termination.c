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

/*************************************************************************
 * IS_PRIMAL_INFEASIBLE
 *************************************************************************/
c_int is_primal_infeasible(QPDOWorkspace *work) {
    size_t n = work->data->n;
    size_t m = work->data->m;
    c_float eps_inf_norm_dy;

    // dy = yh - y
    vec_add_scaled(work->yh, work->y, work->delta_y, -1, m);
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->E, work->delta_y, work->temp_m, m);
        eps_inf_norm_dy = work->settings->eps_prim_inf*vec_norm_inf(work->temp_m, m);
    } else {
        eps_inf_norm_dy = work->settings->eps_prim_inf*vec_norm_inf(work->delta_y, m);
    }

    if (eps_inf_norm_dy == 0) { // dy == 0
        return 0;
    }

    vec_add_scaled(work->Atyh, work->Aty, work->Atdelta_y, -1, n);
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Dinv, work->Atdelta_y, work->Atdelta_y, n);
    }

    // out_of_bounds = u' * max(dy,0) + l' * min(dy,0)
    c_float out_of_bounds = 0;
    if (work->settings->scaling) {
        for (size_t i=0; i < m; i++) {
            out_of_bounds += (work->data->u[i] < work->scaling->E[i]*QPDO_INFTY) ? work->data->u[i]*c_max(work->delta_y[i], 0) : 0;
            out_of_bounds += (work->data->l[i] > -work->scaling->E[i]*QPDO_INFTY) ? work->data->l[i]*c_min(work->delta_y[i], 0) : 0;
        }
    } else {
        for (size_t i=0; i < m; i++) {
            out_of_bounds += (work->data->u[i] < QPDO_INFTY) ? work->data->u[i]*c_max(work->delta_y[i], 0) : 0;
            out_of_bounds += (work->data->l[i] > -QPDO_INFTY) ? work->data->l[i]*c_min(work->delta_y[i], 0) : 0;
        }
    }

    if (vec_norm_inf(work->Atdelta_y, n) <= eps_inf_norm_dy) && (out_of_bounds <= -eps_inf_norm_dy) {
        update_status(work->info, QPDO_PRIMAL_INFEASIBLE);
        if (work->settings->scaling) {
            vec_self_mult_scalar(work->dy, work->scaling->cinv, m);
            vec_ew_prod(work->scaling->E, work->dy, work->dy, m);
        }
        return 1;
    } else {
        return 0;
    }
}

/*************************************************************************
 * IS_DUAL_INFEASIBLE
 *************************************************************************/
c_int is_dual_infeasible(QPDOWorkspace *work) {
    c_float eps_inf_norm_dx, dxQdx, dxdx;
    size_t n = work->data->n;
    size_t m = work->data->m;

    //dx = x-x_prev
    vec_add_scaled(work->x, work->x_prev, work->delta_x, -1, n);
    if (work->settings->scaling) {
        //D*dx
        vec_ew_prod(work->scaling->D, work->delta_x, work->temp_n, n);
        eps_inf_norm_dx = work->settings->eps_dual_inf*vec_norm_inf(work->temp_n, n);
        dxdx = vec_prod(work->temp_n, work->temp_n, n);
    } else {
        eps_inf_norm_dx = work->settings->eps_dual_inf*vec_norm_inf(work->delta_x, n);
        dxdx = vec_prod(work->delta_x, work->delta_x, n);
    }

    if (eps_inf_norm_dx == 0) { //dx == 0
        return 0;
    }

    size_t k;
    //NB Adx = work->Ad (= tau*Ad of the previous iteration)
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Einv, work->Ad, work->Adelta_x, m);
        for (k = 0; k < m; k++) {
            if ((work->data->u[k] < work->scaling->E[k]*QPDO_INFTY && work->Adelta_x[k] >= eps_inf_norm_dx)
                || (work->data->l[k] > -work->scaling->E[k]*QPDO_INFTY && work->Adelta_x[k] <= -eps_inf_norm_dx)) {
                return 0;
            }
        }
    } else {
        for (k = 0; k < m; k++) {
            if ((work->data->u[k] < QPDO_INFTY && work->Ad[k] >= eps_inf_norm_dx)
                || (work->data->l[k] > -QPDO_INFTY && work->Ad[k] <= -eps_inf_norm_dx)) {
                return 0;
            }
        }
    }
    //NB Qdx = work->Qd (= tau*Qd of the previous iteration)
    //NB Qdx = work->Qd - sigma*tau*d (= tau*Qd of the previous iteration) if proximal is used
    if (work->settings->proximal) {
        vec_add_scaled(work->Qdx, work->dx, work->temp_n, - work->sigma * work->tau, n);
        dxQdx = vec_prod(work->delta_x, work->temp_n, n);
    } else {
        dxQdx = vec_prod(work->Qd, work->delta_x, n);
    }
    if (work->settings->scaling) {
        if ((dxQdx <= work->scaling->c*work->settings->eps_dual_inf*work->settings->eps_dual_inf*dxdx)
                && (vec_prod(work->data->q, work->delta_x, n) <= -work->scaling->c*eps_inf_norm_dx)) {
                    update_status(work->info, QPDO_DUAL_INFEASIBLE);
                    vec_ew_prod(work->scaling->D, work->dx, work->dx, n);
                    return 1;
                }
    } else {
        if ((dxQdx <= work->settings->eps_dual_inf*work->settings->eps_dual_inf*dxdx)
                && (vec_prod(work->data->q, work->delta_x, n) <= -eps_inf_norm_dx)) {
                    update_status(work->info, QPDO_DUAL_INFEASIBLE);
                    return 1;
                }
    }
    return 0;
}
