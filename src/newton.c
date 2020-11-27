# ifdef __cplusplus
extern "C" {
# endif

#include "newton.h"
#include "lin_alg.h"
#include "cholmod.h"
#include <stdio.h>

/*************************************************************************
 * NEWTON_DIRECTION
 *************************************************************************/
void newton_direction(QPDOWorkspace *work, cholmod_common *c) {
    // dimensions
    size_t n = work->data->n;
    size_t m = work->data->m;
    // active set
    active_constraints(work);
    // factorization (update)
    enter_leave_constraints(work);
    if ((work->chol->reset_newton && work->chol->n_active_set) || (work->chol->n_enter + work->chol->n_leave) > MAX_RANK_UPDATE) {
        work->chol->reset_newton = FALSE;
        ldlcholQAtmuA(work, c);
    } else if (work->chol->n_active_set) {
        if(work->chol->n_enter) {
            ldlupdate_enter_constraints(work, c);
        }
        if(work->chol->n_leave) {
            ldldowndate_leave_constraints(work, c);
        }
    } else {
        ldlchol(work->data->Q, work, c);
    }
    // right-hand side
    // "dy" = (I+P) .* res_prim_in ./ mu
    // needed later for computing dual search direction dy
    for (size_t i = 0; i < m; i++){
        work->dy[i] = work->res_prim_in[i] / work->mu[i];
        if (!work->chol->active_set[i]) work->dy[i] *= 2;
    }
    // linsys_rhs <- - ( dua_res_in + A' * "dy" )
    mat_tpose_vec(work->data->A, work->chol->dy, work->chol->Atdy, c); // A' * "dy"
    for (size_t i = 0; i < n; i++){
        work->linsys_rhs[i] = - work->res_dual_in[i] - work->Atdy[i];
    }

    // primal search direction
    // dx <- linear system
    ldlsolveLD_rhs(work, c);

    // Qdx = Q*dx
    mat_vec(work->data->Q, work->chol->dx, work->chol->Qdx, c);
    if (work->settings->proximal) {
        vec_add_scaled(work->Qdx, work->dx, work->Qdx, work->sigma, n);
    }
    // Adx = A*dx
    mat_vec(work->data->A, work->chol->dx, work->chol->Adx, c);

    // dual search direction
    // dy <- "dy" + (I-P) .* Adx ./ mu
    for (size_t i = 0; i < m; i++){
       if (work->chol->active_set[i]) work->dy[i] += (work->Adx[i] / work->mu[i]);
    }

    // Atdy = A'*dy
    mat_tpose_vec(work->data->A, work->chol->dy, work->chol->Atdy, c);

    // store old active set
    prea_int_vec_copy(work->chol->active_set, work->chol->active_set_old, m);

    /*//========================= debug mode ============================//
    //=== for checking search direction
    vec_add_scaled(work->Qdx, work->Atdy, work->temp_n, 1, n);
    vec_add_scaled(work->temp_n, work->dua_res_in, work->temp_n, 1, n);
    c_float a1 = vec_norm_inf( work->temp_n, n ); // || (Q + sigma)*dx + A'*dy + dua_res_in ||_inf
    prea_vec_copy(work->pri_res_in, work->temp_m, m);
    size_t i;
    for (i = 0; i < m; i++) {
        if (work->chol->active_set[i]){
            // active : P[i] = 0
            work->temp_m[i] -= work->dy[i] * work->mu[i];
            work->temp_m[i] += work->Adx[i];
        } else {
            // inactive : P[i] = 1
            work->temp_m[i] -= 0.5 * work->dy[i] * work->mu[i];
        }
    }
    c_float a2 = vec_norm_inf( work->temp_m, m ); // || (I-P)*A*dx - (I-P/2)*dy .* mu + pri_res_in ||_inf
    c_print("search direction %.4e %.4e \n", a1, a2);
    //========================= debug mode ============================//*/
}

/*************************************************************************
 * ACTIVE_CONSTRAINTS
 *************************************************************************/
void active_constraints(QPDOWorkspace *work) {
    work->chol->n_active_set = 0;
    for (size_t i = 0; i < work->data->m; i++) {
        if ((work->w[i] <= work->data->l[i]) || ((work->w[i] >= work->data->u[i]))){ // >= and <= ?
            work->chol->active_set[i] = TRUE;
            work->chol->n_active_set++;
        } else {
            work->chol->active_set[i] = FALSE;
        }
    }
    return;
}

/*************************************************************************
 * LEAVE_ENTER_CONSTRAINTS
 *************************************************************************/
void enter_leave_constraints(QPDOWorkspace *work) {
    work->chol->n_enter = 0;
    work->chol->n_leave = 0;
    for (size_t i = 0; i < work->data->m; i++) {
        if (work->chol->active_set[i] && !work->chol->active_set_old[i]) {
            work->chol->enter[ work->chol->n_enter ] = (c_int) i;
            work->chol->n_enter++;
        }
        if (!work->chol->active_set[i] && work->chol->active_set_old[i]) {
            work->chol->leave[ work->chol->n_leave ] = (c_int) i;
            work->chol->n_leave++;
        }
    }
    return;
}

# ifdef __cplusplus
}
# endif
