#include "linesearch.h"
#include "lin_alg.h"
#include <stdlib.h> //for sorting

/*************************************************************************
 * EXACT_LINESEARCH
 *************************************************************************/
void exact_linesearch(QPDOWorkspace *work, cholmod_common *c) {
    //
    // 0.5 * psi^prime( tau ) = eta * tau + beta + delta' * [delta*tau - alpha]_+
    //
    // dimensions
    size_t n = work->data->n;
    size_t m = work->data->m;
    // temp_m = mu .* dy / 2
    vec_ew_prod(work->dy, work->mu, work->temp_m, m);
    vec_self_mult_scalar(work->temp_m, 0.5, m);
    // ls_eta = 0.5 * [ dx' * (Q + sigma * I) * dx + dy' * (mu.*dy)/2 ]
    work->ls_eta = vec_prod(work->dy, work->temp_m, m);
    work->ls_eta += vec_prod(work->dx, work->Qdx, n);
    work->ls_eta *= 0.5;
    // ls_beta = 0.5 * [ dx' * df + y' * (mu.*dy)/2 ]
    // with df = Q*x + q + sigma * (x - xbar)
    work->ls_beta = vec_prod(work->y, work->temp_m, m);
    work->ls_beta += vec_prod(work->dx, work->df, n);
    work->ls_beta *= 0.5;
    // ls_delta = [-c0./sqrt(mu); c0./sqrt(mu)]
    vec_add_scaled(work->Adx, work->temp_m, work->temp_m, -1, m); // c0 = A*dx - (mu.*dy)/2
    vec_ew_prod(work->temp_m, work->sqrt_mu, work->temp_m, m); // c0 ./ sqrt(mu)
    prea_vec_copy(work->temp_m, work->ls_delta + m, m);
    vec_self_mult_scalar(work->temp_m, -1, m);
    prea_vec_copy(work->temp_m, work->ls_delta, m);
    // ls_alpha = [(c1 - l)./sqrt(mu); (u - c1)./sqrt(mu)]
    vec_add_scaled(work->w, work->data->l, work->temp_m, -1, m); // c1 = w, c1 - l
    vec_ew_prod(work->temp_m, work->sqrt_mu, work->temp_m, m); // (c1 - l) ./ sqrt(mu)
    prea_vec_copy(work->temp_m, work->ls_alpha, m);
    vec_add_scaled(work->data->u, work->w, work->temp_m, -1, m); // c1 = w, u - c1
    vec_ew_prod(work->temp_m, work->sqrt_mu, work->temp_m, m); // (u - c1) ./ sqrt(mu)
    prea_vec_copy(work->temp_m, work->ls_alpha + m, m);
    // quit with full step if directional derivative at tau=1 is suff. small
    // // 0.5 * psi^prime( 1 ) = eta + beta + delta' * [delta - alpha]_+
    vec_add_scaled(work->ls_delta, work->ls_alpha, work->temp_2m, -1, m*2);
    size_t i;
    for (i = 0; i < m*2; i++) {
        work->temp_2m[i] = c_max(work->temp_2m[i], 0.0);
    }
    c_float psi_one = vec_prod(work->ls_delta, work->temp_2m, m*2);
    psi_one += work->ls_eta;
    psi_one += work->ls_beta;
    if (c_absval(psi_one) <= 1e-12) {
        return work->tau = 1.0;
    }
    // piecewise-affine linesearch
    pwa_linesearch( work );
    // strictly positive, not crazy stepsize
    //work->tau = c_max( 1e-12, c_min( work->tau, 10.0 ) );
    if (work->tau <= 1e-12) {
        work->tau = 1.0;
    }

    /*//========================= debug mode ============================//
    //=== for checking exact line search
    prea_vec_copy(work->ls_delta, work->temp_2m, m*2);
    vec_self_mult_scalar(work->temp_2m, work->tau, m*2);
    vec_add_scaled(work->temp_2m, work->ls_alpha, work->temp_2m, -1, m*2);
    size_t i;
    for (i = 0; i < m*2; i++) {
        work->temp_2m[i] = c_max( work->temp_2m[i], 0.0 );
    }
    c_float a = vec_prod(work->ls_delta, work->temp_2m, m*2);
    a += work->ls_beta;
    a += work->ls_eta * work->tau;
    c_print("line search %.4e \n", a);
    //========================= debug mode ============================//*/

    return;
}

/*************************************************************************
 * PWA_LINESEARCH
 *************************************************************************/
void pwa_linesearch(QPDOWorkspace *work) {
    //
    // 0 = eta * tau + beta + delta' * [delta*tau - alpha]_+
    //
    // dimensions
    size_t n = work->data->n;
    size_t m = work->data->m;
    // taus = alpha ./ delta
    vec_ew_div(work->ls_alpha, work->ls_delta, work->temp_2m, m*2);
    vec_array_copy(work->temp_2m, work->ls_taus, m*2);
    // idx_L = ( taus > 0 )
    size_t nL = 0;
    size_t i;
    for (i = 0; i < m*2; i++){
        if (work->temp_2m[i] > 0) {
            work->ls_idx_L[i] = TRUE;
            nL++;
        } else {
            work->ls_idx_L[i] = FALSE;
        }
    };

    // taus = taus( idx_L )
    select_subsequence(work->ls_taus, work->ls_taus, work->ls_idx_L, m*2);

    // idx_P = ( delta > 0 )
    for (i = 0; i < m*2; i++){
        if (work->ls_delta[i] > 0) {
            work->ls_idx_P[i] = TRUE;
        } else {
            work->ls_idx_P[i] = FALSE;
        }
    };
    // idx_J = (idx_P & ~idx_L) | (~idx_P & idx_L)
    for (i = 0; i < m*2; i++){
        if ((work->ls_idx_P[i] + work->ls_idx_L[i]) == 1) {
            work->ls_idx_J[i] = TRUE;
        } else {
            work->ls_idx_J[i] = FALSE;
        }
    };

    // a = eta  + delta(idx_J)' * delta(idx_J);
    // b = beta - delta(idx_J)' * alpha(idx_J);
    c_float a, b;
    a = work->ls_eta  + vec_prod_ind(work->ls_delta, work->ls_delta, work->ls_idx_J, m*2);
    b = work->ls_beta - vec_prod_ind(work->ls_delta, work->ls_alpha, work->ls_idx_J, m*2);
    if (nL == 0) {
        work->tau = -b/a;
        return;
    }
    // taus <- sort( taus )
    qsort(work->ls_taus, nL, sizeof(array_element), compare);
    if (b + a * work->ls_taus[0].x > 0) {
        work->tau = -b/a;
        return;
    }
    i = 0;
    size_t iz;
    while (i < nL-1) {
        iz = work->ls_taus[i].i;
        if (work->ls_idx_P[iz]) {
            a = a + work->ls_delta[iz]*work->ls_delta[iz];
            b = b - work->ls_delta[iz]*work->ls_alpha[iz];
        } else {
            a = a - work->ls_delta[iz]*work->ls_delta[iz];
            b = b + work->ls_delta[iz]*work->ls_alpha[iz];
        }
        i++;
        if (b + a * work->ls_taus[i].x > 0) {
            work->tau = -b/a;
            return;
        }
    }
    iz = work->ls_taus[i].i;
    if (work->ls_idx_P[iz]) {
        a = a + work->ls_delta[iz]*work->ls_delta[iz];
        b = b - work->ls_delta[iz]*work->ls_alpha[iz];
    } else {
        a = a - work->ls_delta[iz]*work->ls_delta[iz];
        b = b + work->ls_delta[iz]*work->ls_alpha[iz];
    }
    work->tau = -b/a;
    return;
}

/*************************************************************************
 * VEC_ARRAY_COPY
 *************************************************************************/
void vec_array_copy(c_float *a, array_element* b, size_t n) {
    size_t i;
    array_element ae;
    for (i = 0; i < n; i++) {
        ae.x = a[i];
        ae.i = i;
        b[i] = ae;
    }
}

/*************************************************************************
 * SELECT_SUBSEQUENCE
 *************************************************************************/
void select_subsequence(const array_element *a, array_element *b, const c_int *L, size_t n) {
    size_t i;
    size_t nb_elements = 0;
    for (i = 0; i < n; i++) {
        if (L[i]) {
            b[nb_elements] = a[i];
            nb_elements++;
        }
    }
}

/*************************************************************************
 * VEC_PROD_IND
 *************************************************************************/
c_float vec_prod_ind(const c_float *a, const c_float *b, const c_int *L, size_t n) {
    c_float prod = 0.0;
    size_t i;
    for (i = 0; i < n; i++) {
        if (L[i]) {
            prod += a[i] * b[i];
        }
    }
    return prod;
}

/*************************************************************************
 * COMPARE
 *************************************************************************/
int compare (const void * a, const void * b)
{
    c_float f = ((struct array_element*)a)->x;
    c_float s = ((struct array_element*)b)->x;
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}
