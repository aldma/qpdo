#ifndef LIN_ALG_H
# define LIN_ALG_H

# ifdef __cplusplus
extern "C" {
# endif

# include "types.h"
#include "cholmod_interface.h"

// Copy vector into output
c_float* vec_copy(const c_float *a,
                  size_t         n);

// Copy vector into preallocated vector
void prea_vec_copy(const c_float *a,
                   c_float       *b,
                   size_t         n);

// Copy integer vector into preallocated vector
void prea_int_vec_copy(const c_int *a,
                       c_int       *b,
                       size_t       n);

// Fill float vector with a scalar value
void vec_set_scalar(c_float *a,
                    c_float  sc,
                    size_t   n);

// Fill int vector with a scalar value
void vec_set_scalar_int(c_int *a, 
                        c_int  sc, 
                        size_t n);

// Scalar multiplication of a vector, in place
void vec_self_mult_scalar(c_float *a,
                     c_float  sc,
                     size_t   n);

// Scalar multiplication of a vector
void vec_mult_scalar(const c_float *a,
                     c_float  sc,
                     c_float *b,
                     size_t   n);

// Inner product of two vectors
c_float vec_prod(const c_float *a,
                 const c_float *b,
                 size_t         n);

// 2-norm of a vector
c_float vec_norm_two(const c_float *a, size_t n);

// Infinity norm of a vector
c_float vec_norm_inf(const c_float *a,
                     size_t         n);

// Elementwise addition of a scaled vector to another vector
void vec_add_scaled(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_float        sc,
                    size_t         n);

// Elementwise addition of two scaled vectors
void vec_mult_add_scaled(c_float *a,
                         const c_float *b, 
                         c_float sc1, 
                         c_float sc2, 
                         size_t n);

// Elementwise reciprocal of a vector
void vec_ew_recipr(const c_float *a,
                   c_float       *b,
                   size_t         n);

// Elementwise maximum of two vectors
void vec_ew_max_vec(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    size_t         n);

// Elementwise minimum of two vectors
void vec_ew_min_vec(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    size_t         n);

// Elementwise saturation of a vector
void vec_ew_mid_vec(const c_float *a,
                    const c_float *bmin,
                    const c_float *bmax,
                    c_float       *c,
                    size_t         n);

// Elementwise product of two vectors
void vec_ew_prod(const c_float *a,
                 const c_float *b,
                 c_float       *c,
                 size_t         n);

// Elementwise division of two vectors
void vec_ew_div(const c_float *a,
                const c_float *b,
                c_float       *c,
                size_t         n);

// Elementwise square root of a vector
void vec_ew_sqrt(const c_float *a,
                 c_float       *b,
                 size_t         n);

# ifdef __cplusplus
}
# endif

#endif
