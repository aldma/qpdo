#ifndef UTIL_H
# define UTIL_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"
# include "constants.h"

/*************************************************************************
 * COPY_SETTINGS
 *************************************************************************/
QPDOSettings* copy_settings(const QPDOSettings *settings);

/*************************************************************************
 * C_STRCPY
 * Copy string
 *************************************************************************/
void c_strcpy(char       dest[],
              const char source[]);

/*************************************************************************
 * UPDATE_STATUS
 *************************************************************************/
void update_status(QPDOInfo *info,
                   c_int     status_val);



/*************************************************************************
 * PRINTING FUNCTIONS
 *************************************************************************/
#ifdef PRINTING

void print_header(void);

void print_iteration(c_int iter, QPDOWorkspace *work);

void print_final_message(QPDOWorkspace *work);

#endif


/*************************************************************************
 * TIMER FUNCTIONS
 *************************************************************************/

/*! \cond PRIVATE */

# ifdef PROFILING

// Windows
#  ifdef _WIN32

  // Some R packages clash with elements
  // of the windows.h header, so use a
  // slimmer version for conflict avoidance
  # ifdef R_LANG
    #define NOGDI
  # endif

#   include <windows.h>

struct QPDO_TIMER {
  LARGE_INTEGER tic;
  LARGE_INTEGER toc;
  LARGE_INTEGER freq;
};

// Mac
#  elif defined __APPLE__

#   include <mach/mach_time.h>

/* Use MAC OSX  mach_time for timing */
struct QPDO_TIMER {
  uint64_t                  tic;
  uint64_t                  toc;
  mach_timebase_info_data_t tinfo;
};

// Mac
#  elif defined __MACH__

#   include <mach/mach_time.h>

/* Use MAC OSX  mach_time for timing */
struct QPDO_TIMER {
  uint64_t                  tic;
  uint64_t                  toc;
  mach_timebase_info_data_t tinfo;
};

// Linux
#  elif defined __linux__ // ifdef _WIN32

#   include <time.h>
#   include <sys/time.h>

struct QPDO_TIMER {
  struct timespec tic;
  struct timespec toc;
};

#  endif // ifdef _WIN32

/*! \endcond */


/**
 * Start timer.
 */
void qpdo_tic(QPDOTimer *t);

/**
 * Report time in seconds since last call to qpdo_tic.
 */
c_float qpdo_toc(QPDOTimer *t);

# endif


# ifdef __cplusplus
}
# endif

#endif
