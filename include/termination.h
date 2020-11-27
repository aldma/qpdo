#ifndef TERMINATION_H
# define TERMINATION_H

# ifdef __cplusplus
extern "C" {
# endif

#include "types.h"

c_int check_outer_optimality(QPDOWorkspace *work);

c_int check_inner_optimality(QPDOWorkspace *work);

void compute_inner_residuals_norm(QPDOWorkspace *work, cholmod_common *c);

void compute_outer_residuals_norm(QPDOWorkspace *work);

void store_solution(QPDOWorkspace *work);

# ifdef __cplusplus
}
# endif

#endif
