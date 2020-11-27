#ifndef VALIDATE_H
# define VALIDATE_H

# ifdef __cplusplus
extern "C" {
# endif

# include "types.h"

/*************************************************************************
 * VALIDATE_DATA
 *************************************************************************/
c_int validate_data(const QPDOData *data);

/*************************************************************************
 * VALIDATE_SETTINGS
 *************************************************************************/
c_int validate_settings(const QPDOSettings *settings);


# ifdef __cplusplus
}
# endif

#endif
