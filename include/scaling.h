#ifndef SCALING_H
# define SCALING_H

# ifdef __cplusplus
extern "C" {
# endif

# include "types.h"
# include "lin_alg.h"
# include "constants.h"

/*************************************************************************
 * SCALE DATA
 * Ruiz equilibration of the constraint matrix. Scaling of the cost.
 *************************************************************************/
void scale_data(QPDOWorkspace *work);

# ifdef __cplusplus
}
# endif

#endif
