#ifndef ITERATION_H
#define ITERATION_H

#include "types.h"
#include "global_opts.h"

void compute_inner_residuals(QPDOWorkspace* work, cholmod_common *c);

void compute_outer_residuals(QPDOWorkspace* work, cholmod_common *c);

void update_iterate(QPDOWorkspace *work, cholmod_common *c);

void initialize_mu(QPDOWorkspace *work, cholmod_common *c);

void update_mu(QPDOWorkspace* work, cholmod_common *c);

void update_sigma(QPDOWorkspace* work);

c_float compute_objective(QPDOWorkspace *work);

#endif
