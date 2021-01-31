#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifdef __cplusplus
extern "C" {
# endif

/*************************************************************************
 * BOOLEANS
 *************************************************************************/
#define TRUE 1
#define FALSE 0

/*************************************************************************
 * QPDO STATUSES
 *************************************************************************/
# define QPDO_SOLVED (1)                    /**< problem solved to optimality, given the tolerance */
# define QPDO_DUAL_TERMINATED (2)           /**< problem with dual objective higher than the specified bound */
# define QPDO_NON_CVX (-1)                  /**< problem likely non-convex */
# define QPDO_PRIMAL_DUAL_INFEASIBLE (-2)   /**< problem primal-dual infeasible */
# define QPDO_PRIMAL_INFEASIBLE (-3)        /**< problem primal infeasible  */
# define QPDO_DUAL_INFEASIBLE (-4)          /**< problem dual infeasible */
# define QPDO_MAX_ITER_REACHED (-5)         /**< maximum number of iterations reached */
# define QPDO_MAX_TIME_REACHED (-6)         /**< maximum time exceeded */
# define QPDO_UNSOLVED (-10)                /**< problem unsolved, only setup */
# define QPDO_ERROR (-99)                   /**< error occurred */

/*************************************************************************
 * NULL, NOT-A-NUMBER, INFINITY
 *************************************************************************/
# ifndef QPDO_NULL
    #  define QPDO_NULL 0
# endif
# ifndef QPDO_NAN
    #  define QPDO_NAN ((c_float)0x7fc00000UL)
# endif
# ifndef QPDO_INFTY
    #  define QPDO_INFTY ((c_float)1E+20)
# endif

/*************************************************************************
 * SOLVER PARAMETERS AND SETTINGS
 *************************************************************************/
# define MAX_TIME (QPDO_INFTY)      /**< maximum run time */
# define MAX_ITER (1e4)             /**< maximum number of iterations */
# define INNER_MAX_ITER (1e3)       /**< maximum number of inner iterations */
# define EPS_ABS (1e-6)             /**< absolute tolerance */
# define EPS_ABS_IN (1e0)           /**< inner absolute tolerance */
# define RHO (0.1)                  /**< inner tolerance shrink factor */
# define THETA (0.1)                /**< penalty update criterion parameter */
# define DELTA (1e-2)               /**< penalty update factor */
# define MU_MIN (1e-8)              /**< penalty cap */

# define PROXIMAL (TRUE)            /**< use primal regularization? */
# define SIGMA_INIT (1e-1)          /**< initial primal regularization parameter */
# define SIGMA_UPD (1e-1)           /**< primal regularization parameter shrink factor */
# define SIGMA_MIN (1e-7)           /**< minimum primal regularization parameter */

# define SCALING (10)               /**< scaling iterations */
# define MIN_SCALING (1e-9)         /**< minimum scaling value */
# define MAX_SCALING (1e+9)         /**< maximum scaling value */

# define VERBOSE (TRUE)             /**< print infos? */
# define PRINT_INTERVAL (1)         /**< print every .. iterations */

# define RESET_NEWTON_ITER (100)    /**< re-factorize every .. iterations */
# define MAX_RANK_UPDATE 160        /**< re-factorize if update rank is more than .. */

# define EPS_PRIM_INF (1e-6)        /**< primal infeasibility tolerance */
# define EPS_DUAL_INF (1e-6)        /**< dual infeasibility tolerance */

# ifdef __cplusplus
}
# endif

#endif
