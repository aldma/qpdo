#include "util.h"
#include "lin_alg.h"
#include "global_opts.h"
#include "string.h"

/*************************************************************************
 * UTILITY FUNCTIONS
 *************************************************************************/

/*************************************************************************
 * C_STRCPY
 *************************************************************************/
void c_strcpy(char dest[], const char source[]) {
    size_t i;
    for(i = 0; (dest[i] = source[i]) != '\0'; i++);
}

/*************************************************************************
 * COPY_SETTINGS
 *************************************************************************/
QPDOSettings* copy_settings(const QPDOSettings *settings) {
    QPDOSettings *new = c_malloc(sizeof(QPDOSettings));

    // Copy settings
    new->max_iter                 = settings->max_iter;
    new->inner_max_iter           = settings->inner_max_iter;
    new->eps_abs                  = settings->eps_abs;
    new->eps_abs_in               = settings->eps_abs_in;
    new->eps_prim_inf             = settings->eps_prim_inf;
    new->eps_dual_inf             = settings->eps_dual_inf;
    new->rho                      = settings->rho;
    new->theta                    = settings->theta;
    new->delta                    = settings->delta;
    new->mu_min                   = settings->mu_min;
    new->proximal                 = settings->proximal;
    new->sigma_init               = settings->sigma_init;
    new->sigma_upd                = settings->sigma_upd;
    new->sigma_min                = settings->sigma_min;
    new->scaling                  = settings->scaling;
    new->verbose                  = settings->verbose;
    new->print_interval           = settings->print_interval;
    new->reset_newton_iter        = settings->reset_newton_iter;
    new->max_time                 = settings->max_time;
    return new;
}

/*************************************************************************
 * UPDATE_STATUS
 *************************************************************************/
void update_status(QPDOInfo *info, c_int status_val) {
    // Update status value
    info->status_val = status_val;

    // Update status string depending on status val
    switch (status_val)
    {
    case QPDO_SOLVED:
      c_strcpy(info->status, "solved");
      break;
    case QPDO_DUAL_TERMINATED:
      c_strcpy(info->status, "dual terminated");
      break;
    case QPDO_PRIMAL_INFEASIBLE:
      c_strcpy(info->status, "primal infeasible");
      break;
    case QPDO_DUAL_INFEASIBLE:
      c_strcpy(info->status, "dual infeasible");
      break;
    case QPDO_PRIMAL_DUAL_INFEASIBLE:
      c_strcpy(info->status, "primal-dual infeasible");
      break;
    case QPDO_MAX_TIME_REACHED:
      c_strcpy(info->status, "max time exceeded");
      break;
    case QPDO_MAX_ITER_REACHED:
      c_strcpy(info->status, "maximum iterations reached");
      break;
    case QPDO_UNSOLVED:
      c_strcpy(info->status, "unsolved");
      break;
    case QPDO_ERROR:
      c_strcpy(info->status, "error");
      break;
    default:
      c_strcpy(info->status, "unrecognised status value");
      #ifdef PRINTING
        c_eprint("Unrecognised status value %ld", status_val);
      #endif
      break;
    }
}

/*************************************************************************
 * PRINTING FUNCTIONS
 *************************************************************************/

#ifdef PRINTING
/*************************************************************************
 * PRINT_HEADER
 *************************************************************************/
void print_header(void) {
    c_print("============================================================================\n");
    c_print("===                              QPDO  v0.1                              ===\n");
    c_print("============================================================================\n");
    c_print("  iter |  objective     r.prim     r.dual |  r.p. in    r.d. in   stepsize | \n");
    c_print("============================================================================\n");
}

/*************************************************************************
 * PRINT_ITERATION
 *************************************************************************/
void print_iteration(c_int iter, QPDOWorkspace *work) {
    c_print("%6ld | %+-.3e   %.2e   %.2e | %.2e   %.2e   %.2e | \n", iter,
                     work->info->objective, work->info->res_prim_norm,
                     work->info->res_dual_norm, work->info->res_prim_in_norm,
                     work->info->res_dual_in_norm, work->tau );
}

/*************************************************************************
 * PRINT_FINAL_MESSAGE
 *************************************************************************/
void print_final_message(QPDOWorkspace *work) {
    c_print("============================================================================\n");
    size_t characters_box;
    char buf[80];
    switch (work->info->status_val) {
        case QPDO_SOLVED:
            snprintf(buf, 80, "| QPDO finished successfully.                                              |\n");
            break;
        case QPDO_PRIMAL_INFEASIBLE:
            snprintf(buf, 80, "| QPDO detected a primal infeasible problem.                               |\n");
            break;
        case QPDO_DUAL_INFEASIBLE:
            snprintf(buf, 80, "| QPDO detected a dual infeasible problem.                                 |\n");
            break;
        case QPDO_PRIMAL_DUAL_INFEASIBLE:
            snprintf(buf, 80, "| QPDO detected a primal-dual infeasible problem.                          |\n");
            break;
        case QPDO_MAX_ITER_REACHED:
            snprintf(buf, 80, "| QPDO hit the maximum number of iterations.                               |\n");
            break;
        case QPDO_MAX_TIME_REACHED:
            snprintf(buf, 80, "| QPDO exceeded the specified time limit.                                  |\n");
            break;
        default:
            c_strcpy(work->info->status, "unrecognised status value");
            c_eprint("Unrecognised final status value %ld", work->info->status_val);
            return;
    }
    characters_box = strlen(buf);
    c_print("%s", buf);
    c_print("| primal residual: %5.4e,                primal tolerance: %5.4e |\n", work->info->res_prim_norm, work->settings->eps_abs);
    c_print("| dual residual  : %5.4e,                dual tolerance  : %5.4e |\n", work->info->res_dual_norm, work->settings->eps_abs);
    c_print("| objective value: %+-5.4e                                             |\n", work->info->objective);
    #ifdef PROFILING
        size_t characters_runtime;
        if (work->info->run_time > 1.0) {
            snprintf(buf, 80,"| runtime:         %4.2f seconds", work->info->run_time);
            characters_runtime = strlen(buf);
            c_print("%s", buf);
        } else {
            snprintf(buf, 80,"| runtime:         %4.2f milliseconds", work->info->run_time * 1000);
            characters_runtime = strlen(buf);
            c_print("%s", buf);
        }
        for (; characters_runtime < characters_box-2; characters_runtime++) {
            c_print(" ");
        }
        c_print("|\n");
    #endif
    c_print("============================================================================\n");
    c_print("\n");
}

#endif

/*************************************************************************
 * PROFILING FUNCTIONS
 *************************************************************************/

#ifdef PROFILING

/*************************************************************************
 * QPDO_TIC
 * QPDO_TOC
 *************************************************************************/
// Windows
# ifdef _WIN32
    // tic
    void qpdo_tic(QPDOTimer *t)
    {
        QueryPerformanceFrequency(&t->freq);
        QueryPerformanceCounter(&t->tic);
    }
    // toc
    c_float qpdo_toc(QPDOTimer *t)
    {
        QueryPerformanceCounter(&t->toc);
        return (t->toc.QuadPart - t->tic.QuadPart) / (c_float)t->freq.QuadPart;
    }

// Mac
# elif defined __APPLE__
    // tic
    void qpdo_tic(QPDOTimer *t)
    {
        t->tic = mach_absolute_time();
    }
    // toc
    c_float qpdo_toc(QPDOTimer *t)
    {
        uint64_t duration; // elapsed clock cycles
        t->toc   = mach_absolute_time();
        duration = t->toc - t->tic;
        // from clock cycles to nanoseconds
        mach_timebase_info(&(t->tinfo));
        duration *= t->tinfo.numer;
        duration /= t->tinfo.denom;
        return (c_float)duration / 1e9;
    }

// Mac
# elif defined __MACH__
    // tic
    void qpdo_tic(QPDOTimer *t)
    {
        t->tic = mach_absolute_time();
    }
    // toc
    c_float qpdo_toc(QPDOTimer *t)
    {
        uint64_t duration; // elapsed clock cycles
        t->toc   = mach_absolute_time();
        duration = t->toc - t->tic;
        // from clock cycles to nanoseconds
        mach_timebase_info(&(t->tinfo));
        duration *= t->tinfo.numer;
        duration /= t->tinfo.denom;
        return (c_float)duration / 1e9;
    }

// Linux
# elif defined __linux__
    // tic
    void qpdo_tic(QPDOTimer *t)
    {
        clock_gettime(CLOCK_MONOTONIC, &t->tic);
    }
    // toc
    c_float qpdo_toc(QPDOTimer *t)
    {
        struct timespec temp;

        clock_gettime(CLOCK_MONOTONIC, &t->toc);

        if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
            temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec - 1;
            temp.tv_nsec = 1000000000 + t->toc.tv_nsec - t->tic.tv_nsec;
        } else {
            temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec;
            temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
        }
        return (c_float)temp.tv_sec + (c_float)temp.tv_nsec / 1e9;
    }

# endif

#endif
