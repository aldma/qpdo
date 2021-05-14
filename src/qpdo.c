# ifdef __cplusplus
extern "C" {
# endif

#include "qpdo.h"
#include "global_opts.h"
#include "constants.h"
#include "validate.h"
#include "lin_alg.h"
#include "util.h"
#include "scaling.h"
#include "linesearch.h"
#include "termination.h"
#include "cholmod.h"
#include "cholmod_function.h"
#include "cholmod_interface.h"
#include "newton.h"
#include "iteration.h"

/*************************************************************************
 * QPDO_SET_DEFAULT_SETTINGS
 * refer to constants.h
 *************************************************************************/
void qpdo_set_default_settings(QPDOSettings *settings) {
    settings->max_time                = MAX_TIME;                /* maximum run time */
    settings->max_iter                = MAX_ITER;                /* maximum iterations */
    settings->inner_max_iter          = INNER_MAX_ITER;          /* maximum inner iterations */
    settings->eps_abs                 = (c_float)EPS_ABS;        /* outer absolute tolerance */
    settings->eps_abs_in              = (c_float)EPS_ABS_IN;     /* inner absolute tolerance */
    settings->eps_prim_inf            = (c_float)EPS_PRIM_INF;   /* primal infeasibility tolerance */
    settings->eps_dual_inf            = (c_float)EPS_DUAL_INF;   /* dual infeasibility tolerance */
    settings->rho                     = (c_float)RHO;            /* inner tolerance shrink factor */
    settings->theta                   = (c_float)THETA;          /* penalty update check factor */
    settings->delta                   = (c_float)DELTA;          /* penalty update shrink factor */
    settings->mu_min                  = (c_float)MU_MIN;         /* minimum dual regularization/penalty parameter */
    settings->proximal                = PROXIMAL;                /* boolean, use primal regularization? */
    settings->sigma_init              = (c_float)SIGMA_INIT;     /* initial primal regularization parameter */
    settings->sigma_upd               = (c_float)SIGMA_UPD;      /* primal regularization parameter shrink factor */
    settings->sigma_min               = (c_float)SIGMA_MIN;      /* minimum primal regularization parameter */
    settings->scaling                 = SCALING;                 /* int, scaling */
    settings->verbose                 = VERBOSE;                 /* boolean, print infos? */
    settings->print_interval          = PRINT_INTERVAL;          /* int, print every .. iterations */
    settings->reset_newton_iter       = RESET_NEWTON_ITER;       /* int, re-factorize every .. iterations */
}

/*************************************************************************
 * QPDO_SETUP
 *************************************************************************/
QPDOWorkspace* qpdo_setup(const QPDOData *data, const QPDOSettings *settings) {

    QPDOWorkspace *work;

    // Validate data
    if (!validate_data(data)) {
        # ifdef PRINTING
        c_eprint("Data validation returned failure");
        # endif
        return QPDO_NULL;
    }

  // Validate settings
    if (!validate_settings(settings)) {
        # ifdef PRINTING
        c_eprint("Settings validation returned failure");
        # endif
        return QPDO_NULL;
    }

    // Allocate empty workspace
    work = c_calloc(1, sizeof(QPDOWorkspace));
    if (!work) {
        # ifdef PRINTING
        c_eprint("allocating work failure");
        # endif
        return QPDO_NULL;
    }

    // Start and allocate directly timer
    # ifdef PROFILING
    work->timer = c_malloc(sizeof(QPDOTimer));
    qpdo_tic(work->timer);
    # endif

    // Copy settings
    work->settings = copy_settings(settings);
    work->sqrt_delta = c_sqrt(work->settings->delta);
    work->sigma = work->settings->sigma_init;

    size_t n = data->n;
    size_t m = data->m;

    // initialize CHOLMOD and its settings
    work->chol = c_calloc(1, sizeof(QPDOCholmod));
    cholmod_common common, *c;
    c = &common;
    CHOLMOD(start)(c);
    cholmod_set_settings(c);

    // Copy problem data into workspace
    work->data      = c_calloc(1, sizeof(QPDOData));
    work->data->n   = data->n;
    work->data->m   = data->m;
    work->data->l   = vec_copy(data->l, m);
    work->data->u   = vec_copy(data->u, m);
    work->data->q   = vec_copy(data->q, n);
    work->data->c   = data->c;
    work->data->Q   = CHOLMOD(copy_sparse)(data->Q, c);
    work->data->A   = CHOLMOD(copy_sparse)(data->A, c);
    work->data->A->stype = 0;

    // Allocate internal solver variables
    work->x     = c_calloc(n, sizeof(c_float));
    work->y     = c_calloc(m, sizeof(c_float));
    work->xbar  = c_calloc(n, sizeof(c_float));
    work->ybar  = c_calloc(m, sizeof(c_float));
    work->Ax    = c_calloc(m, sizeof(c_float));
    work->Qx    = c_calloc(n, sizeof(c_float));
    work->Aty   = c_calloc(n, sizeof(c_float));

    // Workspace variables
    work->initialized       = FALSE;
    work->temp_m            = c_calloc(m, sizeof(c_float));
    work->temp_n            = c_calloc(n, sizeof(c_float));
    work->temp_2m           = c_calloc(m*2, sizeof(c_float));
    work->mu                = c_calloc(m, sizeof(c_float));
    work->sqrt_mu           = c_calloc(m, sizeof(c_float));
    work->n_mu_changed      = 0;

    // compute_residuals
    work->z             = c_calloc(m, sizeof(c_float));
    work->w             = c_calloc(m, sizeof(c_float));
    work->res_prim      = c_calloc(m, sizeof(c_float));
    work->res_prim_old  = c_calloc(m, sizeof(c_float));
    work->res_dual      = c_calloc(n, sizeof(c_float));
    work->res_prim_in   = c_calloc(m, sizeof(c_float));
    work->res_dual_in   = c_calloc(n, sizeof(c_float));
    work->df            = c_calloc(n, sizeof(c_float));

    // Linesearch variables
    work->ls_delta    = c_calloc(m*2, sizeof(c_float));
    work->ls_alpha    = c_calloc(m*2, sizeof(c_float));
    work->ls_taus     = c_calloc(m*2, sizeof(array_element));
    work->ls_idx_L    = c_calloc(m*2, sizeof(c_int));
    work->ls_idx_P    = c_calloc(m*2, sizeof(c_int));
    work->ls_idx_J    = c_calloc(m*2, sizeof(c_int));

    // Perform scaling
    if (settings->scaling) {
        // Allocate scaling structure
        work->scaling       = c_malloc(sizeof(QPDOScaling));
        work->scaling->D    = c_calloc(n, sizeof(c_float));
        work->scaling->Dinv = c_calloc(n, sizeof(c_float));
        work->scaling->E    = c_calloc(m, sizeof(c_float));
        work->scaling->Einv = c_calloc(m, sizeof(c_float));
        // Allocate cholmod_dense pointers to E_temp and D_temp
        work->chol->E_temp = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, c);
        work->E_temp = work->chol->E_temp->x;
        work->chol->D_temp = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
        work->D_temp = work->chol->D_temp->x;
        // Scale data
        scale_data(work);
        // || q ||
        vec_ew_prod(work->scaling->Dinv, work->data->q, work->temp_n, n);
        work->norm_q = vec_norm_inf(work->temp_n, n);
    }
    else {
        work->scaling = QPDO_NULL;
        // || q ||
        work->norm_q = vec_norm_inf(work->data->q, n);
    }

    // CHOLMOD variables
    work->chol->linsys_rhs = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
    work->linsys_rhs    = work->chol->linsys_rhs->x;
    work->chol->dx      = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
    work->dx            = work->chol->dx->x;
    work->chol->dy      = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, c);
    work->dy            = work->chol->dy->x;
    work->chol->Qdx     = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
    work->Qdx           = work->chol->Qdx->x;
    work->chol->Adx     = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, c);
    work->Adx           = work->chol->Adx->x;
    work->chol->Atdy    = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
    work->Atdy          = work->chol->Atdy->x;
    work->chol->active_set      = c_calloc(m, sizeof(c_int));
    work->chol->active_set_old  = c_calloc(m, sizeof(c_int));
    vec_set_scalar_int(work->chol->active_set_old, FALSE, m);
    work->chol->reset_newton = TRUE;
    work->chol->enter       = c_calloc(m, sizeof(c_int));
    work->chol->leave       = c_calloc(m, sizeof(c_int));
    work->chol->At_scale    = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, c);

    // Allocate solution
    work->solution    = c_calloc(1, sizeof(QPDOSolution));
    work->solution->x = c_calloc(1, n * sizeof(c_float));
    work->solution->y = c_calloc(1, m * sizeof(c_float));

    // Allocate and initialize information
    work->info        = c_calloc(1, sizeof(QPDOInfo));
    update_status(work->info, QPDO_UNSOLVED);
    # ifdef PROFILING
    work->info->solve_time  = 0.0;                   // Solve time to zero
    work->info->run_time    = 0.0;                   // Total run time to zero
    work->info->setup_time  = qpdo_toc(work->timer); // Update timer information
    # endif

    //Finish cholmod
    CHOLMOD(finish)(c);

    // Return workspace structure
    return work;
}

/*************************************************************************
 * QPDO_WARM_START
 *************************************************************************/
void qpdo_warm_start(QPDOWorkspace *work, c_float *x_warm_start, c_float *y_warm_start) {

    // initialize primal regularization
    work->sigma = work->settings->sigma_init;

    // if problem not UNSOLVED anymore, then just count the warm start as the setup time
    #ifdef PROFILING
    if (work->info->status_val != QPDO_UNSOLVED) work->info->setup_time = 0;
    qpdo_tic(work->timer);
    #endif

    size_t n = work->data->n;
    size_t m = work->data->m;

    cholmod_common common, *c;
    c = &common;
    CHOLMOD(start)(c);
    cholmod_set_settings(c);

    // warm-start primal variable
    if (x_warm_start != NULL) {
        prea_vec_copy(x_warm_start, work->x, n);
        // scaling
        if (work->settings->scaling) {
            vec_ew_prod(work->x, work->scaling->Dinv, work->x, n);
        }
        prea_vec_copy(work->x, work->xbar, n);

        // NB: link to Cholmod using "dx" and "chol->dx" as "x"
        prea_vec_copy(work->x, work->dx, n);
        mat_vec(work->data->Q, work->chol->dx, work->chol->Qdx, c);
        if (work->settings->proximal) {
            vec_add_scaled(work->Qdx, work->x, work->Qx, work->sigma, n);
        } else {
            prea_vec_copy(work->Qdx, work->Qx, n);
        }
        mat_vec(work->data->A, work->chol->dx, work->chol->Adx, c);
        prea_vec_copy(work->Adx, work->Ax, m);

        work->info->objective = compute_objective(work);

    } else {
        // if not warm-started, set equal to zero
        vec_set_scalar(work->x, 0., n);
        vec_set_scalar(work->xbar, 0., n);
        vec_set_scalar(work->Qx, 0., n);
        vec_set_scalar(work->Ax, 0., m);
        work->info->objective = 0.0;
    }

    // warm-start dual variable
    if (y_warm_start != NULL) {
        prea_vec_copy(y_warm_start, work->y, m);
        // scaling
        if (work->settings->scaling) {
            vec_ew_prod(work->y, work->scaling->Einv, work->y, m);
            vec_self_mult_scalar(work->y, work->scaling->c, m);
        }
        prea_vec_copy(work->y, work->ybar, m);

        prea_vec_copy(work->y, work->dy, m);
        mat_tpose_vec(work->data->A, work->chol->dy, work->chol->Atdy, c);
        prea_vec_copy(work->Atdy, work->Aty, n);

    } else {
        // if not warm-started, set equal to zero
        vec_set_scalar(work->y, 0., m);
        vec_set_scalar(work->ybar, 0., m);
        vec_set_scalar(work->Aty, 0., n);
    }

    // initialize dual regularization
    initialize_mu(work, c);

    work->initialized = TRUE;

    CHOLMOD(finish)(c);

    // update timer
    #ifdef PROFILING
    work->info->setup_time += qpdo_toc(work->timer);
    #endif
}

/*************************************************************************
 * QPDO_SOLVE
 *************************************************************************/
void qpdo_solve(QPDOWorkspace *work) {

    // print header
    #ifdef PRINTING
    if (work->settings->verbose) print_header();
    #endif

    // check if workspace was correctly initialized
    if (!work->initialized) {
        qpdo_warm_start(work, NULL, NULL);
    }

    // dimensions
    size_t n = work->data->n;
    size_t m = work->data->m;

    //Set initial workspace variables
    work->eps_in = work->settings->eps_abs_in;
    work->sigma = work->settings->sigma_init;
    work->chol->reset_newton = TRUE;
    vec_set_scalar_int(work->chol->active_set_old, FALSE, m);

    // Start the timer (after warm_start because this is already added to the setup time)
    #ifdef PROFILING
    qpdo_tic(work->timer);
    #endif

    // initialize CHOLMOD and its settings
    cholmod_common common, *c;
    c = &common;
    CHOLMOD(start)(c);
    cholmod_set_settings(c);

    // counters
    c_int iter; // iteration counter
    c_int oter = 0; // outer iteration counter (number of subproblems)
    c_int iter_old = 0; // iteration number ´iter´ at which the previous subproblem finished

    // main loop
    for (iter = 0; iter < work->settings->max_iter; iter++) {

        // compute outer residuals and norm
        compute_outer_residuals(work, c);
        compute_outer_residuals_norm(work);

        // compute inner residuals and norm
        compute_inner_residuals(work, c);
        compute_inner_residuals_norm( work, c);

        #ifdef PRINTING
        if ((work->settings->verbose) && (iter % work->settings->print_interval == 0)) {
            work->info->objective = compute_objective(work);
            print_iteration(iter, work);
        }
        #endif

        // check optimality
        if (check_outer_optimality(work)) {
            break; // problem solved
        }

        if (((iter > iter_old + 1) && check_inner_optimality(work)) || (iter == iter_old + work->settings->inner_max_iter)) {
            // subproblem solved or stopped

            if (iter < iter_old + work->settings->inner_max_iter) {

                if (work->settings->eps_prim_inf > 0) {
                    // dy = y - ybar
                    vec_add_scaled(work->y, work->ybar, work->dy, -1, m);
                    // Atdy = A' * dy
                    mat_tpose_vec(work->data->A, work->chol->dy, work->chol->Atdy, c);
                    // check primal infeasibility
                    if (is_primal_infeasible(work)) {
                        break;
                    }
                }

                if (work->settings->eps_dual_inf > 0) {
                    // dx = x - xbar
                    vec_add_scaled(work->x, work->xbar, work->dx, -1, n);
                    // Qdx = Q * dx
                    mat_vec(work->data->Q, work->chol->dx, work->chol->Qdx, c);
                    // Adx = A * dx
                    mat_vec(work->data->A, work->chol->dx, work->chol->Adx, c);
                    // check dual infeasibility
                    if (is_dual_infeasible(work)) {
                        break;
                    }
                }
            }

            // update iterate estimate
            prea_vec_copy(work->x, work->xbar, n);
            prea_vec_copy(work->y, work->ybar, m);

            // update regularization parameters
            if ((oter > 0) && (work->info->res_prim_norm > work->settings->eps_abs)) {
                update_mu(work, c);
            }
            if ((work->settings->proximal) && (oter > 0) && (work->info->res_dual_norm > work->settings->eps_abs)) {
                update_sigma(work);
            }

            if (iter < iter_old + work->settings->inner_max_iter) {

                // update inner tolerance
                work->eps_in = c_max( work->settings->rho * work->eps_in , work->settings->eps_abs );

                #ifdef PRINTING
                if ((work->settings->verbose) && (iter % work->settings->print_interval == 0)) {
                    c_print("%6ld |-------------------------------------------------------------------|\n", iter);
                }
                #endif
            } else {
                #ifdef PRINTING
                if ((work->settings->verbose) && (iter % work->settings->print_interval == 0)) {
                    c_print("%6ld |--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- |\n", iter);
                }
                #endif
            }

            // shift old primal residual
            prea_vec_copy(work->res_prim, work->res_prim_old, m);

            // update counters
            oter++;
            iter_old = iter;

        } else {
            // subproblem iteration
            if (iter % work->settings->reset_newton_iter == 0) {
                work->chol->reset_newton = TRUE;
            }
            // primal-dual step
            update_iterate(work, c);
        }

        // check max time reached
        #ifdef PROFILING
        work->info->run_time = work->info->setup_time + qpdo_toc(work->timer);
        if (work->info->run_time > work->settings->max_time) {
            update_status(work->info, QPDO_MAX_TIME_REACHED);
            break; // max time reached
        }
        #endif
    }

    if (work->info->status_val == QPDO_UNSOLVED) {
        update_status(work->info, QPDO_MAX_ITER_REACHED); // max iterations reached
    }

    work->info->iterations = iter;
    work->info->oterations = oter;
    store_solution(work);
    CHOLMOD(finish)(c);
    work->initialized = FALSE;
    // profiling
    #ifdef PROFILING
    work->info->solve_time = qpdo_toc(work->timer);
    work->info->run_time = work->info->setup_time + work->info->solve_time;
    #endif
    // print final message
    #ifdef PRINTING
    if (work->settings->verbose) {
        if (iter % work->settings->print_interval != 0) {
            print_iteration(iter, work);
        }
        print_final_message(work);
    }
    #endif

    return;
} // qpdo_solve

/*************************************************************************
 * QPDO_UPDATE_SETTINGS
 *************************************************************************/
void qpdo_update_settings(QPDOWorkspace* work, const QPDOSettings *settings) {
    // Validate settings
    if (!validate_settings(settings)) {
        # ifdef PRINTING
        c_eprint("Settings validation returned failure");
        # endif
        update_status(work->info, QPDO_ERROR);
        return;
    }
    if (work->settings->scaling > settings->scaling) {
        # ifdef PRINTING
        c_eprint("Decreasing the number of scaling iterations is not allowed");
        # endif
        update_status(work->info, QPDO_ERROR);
        return;
    } else if (work->settings->scaling < settings->scaling) {
        // Save current scaling vectors
        prea_vec_copy(work->scaling->D, work->temp_n, work->data->n);
        prea_vec_copy(work->scaling->E, work->temp_m, work->data->m);
        c_float c_temp = work->scaling->c;
        // Perform the remaining scaling iterations
        work->settings->scaling = settings->scaling - work->settings->scaling;
        scale_data(work);
        // Compute the total scaling vectors
        vec_ew_prod(work->scaling->D, work->temp_n, work->scaling->D, work->data->n);
        vec_ew_prod(work->scaling->E, work->temp_m, work->scaling->E, work->data->m);
        work->scaling->c *= c_temp;
        // Save the inverses
        vec_ew_recipr(work->scaling->D, work->scaling->Dinv, work->data->n);
        vec_ew_recipr(work->scaling->E, work->scaling->Einv, work->data->m);
        work->scaling->cinv = 1/work->scaling->c;
    }
    // Copy settings
    c_free(work->settings);
    work->settings = copy_settings(settings);
    work->sqrt_delta = c_sqrt( work->settings->delta );
}

/*************************************************************************
 * QPDO_UPDATE_BOUNDS
 *************************************************************************/
void qpdo_update_bounds(QPDOWorkspace *work, const c_float *l, const c_float *u) {
    // Validate bounds
    size_t j;
    size_t m = work->data->m;
    if (l != NULL && u != NULL) {
        for (j = 0; j < m; j++) {
            if (l[j] > u[j]) {
                # ifdef PRINTING
                c_eprint("Lower bound at index %d is greater than upper bound: %.4e > %.4e",
                                        (int)j, work->data->l[j], work->data->u[j]);
                # endif
                update_status(work->info, QPDO_ERROR);
                return;
            }
        }
    }
    if (l != NULL) prea_vec_copy(l, work->data->l, m);
    if (u != NULL) prea_vec_copy(u, work->data->u, m);
    if (work->settings->scaling) {
        if (l != NULL) vec_ew_prod(work->scaling->E, work->data->l, work->data->l, m);
        if (u != NULL) vec_ew_prod(work->scaling->E, work->data->u, work->data->u, m);
    }
}

/*************************************************************************
 * QPDO_UPDATE_Q
 *************************************************************************/
void qpdo_update_q(QPDOWorkspace *work, const c_float *q) {
    size_t n = work->data->n;
    prea_vec_copy(q, work->data->q, n);
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->D, work->data->q, work->data->q, n);
        // Update cost scaling scalar
        c_float c_old = work->scaling->c;
        if (work->settings->proximal) {
            vec_add_scaled(work->Qx, work->x, work->Qx, -work->sigma, n);
        }
        vec_add_scaled(work->data->q, work->Qx, work->temp_n, work->scaling->cinv, n);
        work->scaling->c = 1/c_max(1.0, vec_norm_inf(work->temp_n, n));
        work->scaling->cinv = 1/work->scaling->c;
        vec_self_mult_scalar(work->data->q, work->scaling->c, n);
        // Cholmod
        cholmod_common common, *c;
        c = &common;
        CHOLMOD(start)(c);
        cholmod_dense *scalar = CHOLMOD(ones)(1,1,CHOLMOD_REAL, c);
        c_float *scalarx = scalar->x;
        scalarx[0] = work->scaling->c/c_old;
        CHOLMOD(scale)(scalar, CHOLMOD_SCALAR, work->data->Q, c);
        CHOLMOD(free_dense)(&scalar, c);
        CHOLMOD(finish)(c);
        // adjust Qx
        vec_self_mult_scalar(work->Qx, work->scaling->c/c_old, n);
        if (work->settings->proximal) {
            work->sigma = work->settings->sigma_init;
            vec_add_scaled(work->Qx, work->x, work->Qx, work->sigma, n);
        }
        // adjust || q ||
        vec_ew_prod(work->scaling->Dinv, work->data->q, work->temp_n, n);
        work->norm_q = vec_norm_inf(work->temp_n, n);
    } else {
        // adjust || q ||
        work->norm_q = vec_norm_inf(work->data->q, n);
    }
}

/*************************************************************************
 * QPDO_CLEANUP
 *************************************************************************/
void qpdo_cleanup(QPDOWorkspace *work) {
    if (work) { // If workspace has been allocated

        // Free Data
        cholmod_common common, *c;
        c = &common;
        if (work->data) {
            CHOLMOD(start)(c);
            if (work->data->Q) CHOLMOD(free_sparse)(&work->data->Q, c);
            if (work->data->A) CHOLMOD(free_sparse)(&work->data->A, c);
            CHOLMOD(finish)(c);
            if (work->data->q) c_free(work->data->q);
            if (work->data->l) c_free(work->data->l);
            if (work->data->u) c_free(work->data->u);
            c_free(work->data);
        }

        // Free scaling
        if (work->settings->scaling) {
            if (work->scaling->D) c_free(work->scaling->D);
            if (work->scaling->Dinv) c_free(work->scaling->Dinv);
            if (work->scaling->E) c_free(work->scaling->E);
            if (work->scaling->Einv) c_free(work->scaling->Einv);
            c_free(work->scaling);
        }

        // Free other Variables
        if (work->x) c_free(work->x);
        if (work->y) c_free(work->y);
        if (work->Ax) c_free(work->Ax);
        if (work->Qx) c_free(work->Qx);
        if (work->Aty) c_free(work->Aty);
        if (work->temp_m) c_free(work->temp_m);
        if (work->temp_n) c_free(work->temp_n);
        if (work->temp_2m) c_free(work->temp_2m);
        if (work->mu) c_free(work->mu);
        if (work->z) c_free(work->z);
        if (work->w) c_free(work->w);
        if (work->res_prim) c_free(work->res_prim);
        if (work->res_dual) c_free(work->res_dual);
        if (work->res_prim_old) c_free(work->res_prim_old);
        if (work->res_prim_in) c_free(work->res_prim_in);
        if (work->res_dual_in) c_free(work->res_dual_in);
        if (work->df) c_free(work->df);
        if (work->xbar) c_free(work->xbar);
        if (work->ybar) c_free(work->ybar);
        if (work->sqrt_mu) c_free(work->sqrt_mu);
        if (work->ls_delta) c_free(work->ls_delta);
        if (work->ls_alpha) c_free(work->ls_alpha);
        if (work->ls_taus) c_free(work->ls_taus);
        if (work->ls_idx_L) c_free(work->ls_idx_L);
        if (work->ls_idx_P) c_free(work->ls_idx_P);
        if (work->ls_idx_J) c_free(work->ls_idx_J);

        // Free Settings
        if (work->settings) c_free(work->settings);

        //Free chol struct
        if (work->chol) {
            CHOLMOD(start)(c);
            if (work->chol->D_temp) CHOLMOD(free_dense)(&work->chol->D_temp, c);
            if (work->chol->E_temp) CHOLMOD(free_dense)(&work->chol->E_temp, c);
            if (work->chol->linsys_rhs) CHOLMOD(free_dense)(&work->chol->linsys_rhs, c);
            if (work->chol->dx) CHOLMOD(free_dense)(&work->chol->dx, c);
            if (work->chol->dy) CHOLMOD(free_dense)(&work->chol->dy, c);
            if (work->chol->Qdx) CHOLMOD(free_dense)(&work->chol->Qdx, c);
            if (work->chol->Adx) CHOLMOD(free_dense)(&work->chol->Adx, c);
            if (work->chol->Atdy) CHOLMOD(free_dense)(&work->chol->Atdy, c);
            if (work->chol->LD) CHOLMOD(free_factor)(&work->chol->LD, c);
            if (work->chol->LD_Q) CHOLMOD(free_factor)(&work->chol->LD_Q, c);
            if (work->chol->active_set) c_free(work->chol->active_set);
            if (work->chol->active_set_old) c_free(work->chol->active_set_old);
            if (work->chol->enter) c_free(work->chol->enter);
            if (work->chol->leave) c_free(work->chol->leave);
            if (work->chol->At_scale) CHOLMOD(free_dense)(&work->chol->At_scale, c);
            if (work->chol->At_sqrt_mu) CHOLMOD(free_sparse)(&work->chol->At_sqrt_mu, c);
            CHOLMOD(finish)(c);
            c_free(work->chol);
        }

        // Free solution
        if (work->solution) {
            if (work->solution->x) c_free(work->solution->x);
            if (work->solution->y) c_free(work->solution->y);
            c_free(work->solution);
        }

        // Free timer
        # ifdef PROFILING
            if (work->timer) c_free(work->timer);
        # endif

        // Free information
        if (work->info) c_free(work->info);

        // Free work
        c_free(work);
    }
}

# ifdef __cplusplus
}
# endif
