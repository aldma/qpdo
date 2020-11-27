#include "validate.h"
#include "lin_alg.h"
#include "constants.h"

/***********************************************************
* Validation of data and settings
***********************************************************/

c_int validate_data(const QPDOData *data) {


  if (!data) {
# ifdef PRINTING
    c_eprint("Missing data");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  // Lower and upper bounds
  size_t j;
  for (j = 0; j < data->m; j++) {
    if (data->l[j] > data->u[j]) {
# ifdef PRINTING
      c_eprint("Lower bound at index %d is greater than upper bound: %.4e > %.4e",
               (int)j, data->l[j], data->u[j]);
# endif /* ifdef PRINTING */
      return FALSE;
    }
  }
  return TRUE;
}


c_int validate_settings(const QPDOSettings *settings) {

  if (!settings) {
# ifdef PRINTING
    c_eprint("Missing settings!");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->max_iter <= 0) {
# ifdef PRINTING
    c_eprint("max_iter must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->inner_max_iter <= 0) {
# ifdef PRINTING
    c_eprint("inner_max_iter must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->eps_abs <= 0) {
# ifdef PRINTING
    c_eprint("eps_abs must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

    if (settings->eps_abs_in <= 0) {
# ifdef PRINTING
    c_eprint("eps_abs_in must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->rho <= 0 || settings->rho >= 1) {
# ifdef PRINTING
    c_eprint("rho must be positive and smaller than 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->theta <= 0 || settings->theta > 1) {
# ifdef PRINTING
    c_eprint("theta must be positive and smaller than or equal to 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

    if (settings->delta <= 0 || settings->delta >= 1) {
# ifdef PRINTING
    c_eprint("delta must be positive and smaller than 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->mu_min <= 0) {
# ifdef PRINTING
    c_eprint("mu_min must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if ((settings->proximal != 0) && (settings->proximal != 1)) {
# ifdef PRINTING
    c_eprint("proximal must be either 0 or 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if (settings->sigma_init <= 0) {
# ifdef PRINTING
    c_eprint("sigma_init must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if (settings->sigma_upd <= 0 || settings->sigma_upd > 1) {
# ifdef PRINTING
    c_eprint("sigma_upd must be positive and smaller than or equal to 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if (settings->sigma_min > settings->sigma_init) {
# ifdef PRINTING
    c_eprint("sigma_min must be smaller than or equal to sigma_init");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if (settings->scaling < 0) {
# ifdef PRINTING
    c_eprint("scaling must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->verbose < 0) {
# ifdef PRINTING
    c_eprint("verbose must be either nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->print_interval < 0) {
# ifdef PRINTING
    c_eprint("print_interval must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->reset_newton_iter < 0) {
# ifdef PRINTING
    c_eprint("reset_newton_iter must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  return TRUE;
}
