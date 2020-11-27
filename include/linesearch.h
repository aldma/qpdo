#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "types.h"
#include "constants.h"

// Perform exact linesearch
void exact_linesearch(QPDOWorkspace *work, cholmod_common *c);

// Root-finder of piecewise affine function (using qsort)
void pwa_linesearch(QPDOWorkspace *work);

// Helper function to copy vector into array
void vec_array_copy(c_float       *a,
                    array_element *b,
                    size_t         n);

// Select elements of array based on a set of indices
void select_subsequence(const array_element *a,
                        array_element       *b,
                        const c_int         *L,
                        size_t               n);

// Inner product over index set
c_float vec_prod_ind(const c_float *a,
                     const c_float *b,
                     const c_int   *L,
                     size_t         n);

// Helper function for qsort
int compare (const void * a,
             const void * b);

#endif
