#include "qpdo.h"
#include <stdio.h>

#define N 2
#define M 3
#define ANZMAX 4
#define QNZMAX 2

c_float* random_vector(c_int n) {
  c_float* X = c_calloc(n, sizeof(c_float));
  for (int i = 0; i < n; i++) {
    X[i] = (c_float) 10*rand()/RAND_MAX;
  }
  return X;
}

c_float* constant_vector(c_float c, c_int n) {
  c_float* X = c_calloc(n, sizeof(c_float));
  for (int i = 0; i < n; i++) {
    X[i] = c;
  }
  return X;
}

int main() {

    // Problem data
    QPDOData *data;
    data    = (QPDOData *)c_malloc(sizeof(QPDOData));
    data->n = (c_int) N;
    data->m = (c_int) M;
    data->q = random_vector(data->n);
    data->c = 0;
    data->l = constant_vector(-2, data->m);
    data->u = constant_vector(2, data->m);

    // Workspace
    QPDOWorkspace *work;

    // Settings
    QPDOSettings *settings = (QPDOSettings *)c_malloc(sizeof(QPDOSettings));

    // Linear solver setup
    solver_common common, *c;
    c = &common;
    solver_sparse *A, *Q;

    #include "cholmod.h"
    CHOLMOD(start)(c);
    A = CHOLMOD(allocate_sparse)(m, n, ANZMAX, 1, 1, 0, CHOLMOD_REAL, c);
    Q = CHOLMOD(allocate_sparse)(n, n, QNZMAX, 1, 1, -1, CHOLMOD_REAL, c);

    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 1.0; Ax[1] = 1.0; Ax[2] = 1.0; Ax[3] = 1.0;
    Ap[0] = 0; Ap[1] = 2; Ap[2] = 4;
    Ai[0] = 0; Ai[1] = 2; Ai[2] = 1; Ai[3] = 2;
    c_float *Qx;
    c_int *Qi, *Qp;
    Qx = Q->x;
    Qp = Q->p;
    Qi = Q->i;
    Qx[0] = 1.0; Qx[1] = 1.5;
    Qp[0] = 0; Qp[1] = 1; Qp[2] = 2;
    Qi[0] = 0; Qi[1] = 1;

    data->A = A;
    data->Q = Q;

    // Set default settings
    qpdo_set_default_settings(settings);

    // Workspace setup
    work = qpdo_setup(data, settings, c);

    // Solve
    qpdo_solve(work);

    printf("Solver status: ");
    printf(work->info->status);
    printf(" \n");
    printf("Iterations: %d\n", (int) work->info->iterations);
    printf("Outer iter: %d\n", (int) work->info->oterations);

    // Clean and free memory
    CHOLMOD(start)(c);
    CHOLMOD(free_sparse)(&data->Q, c);
    CHOLMOD(free_sparse)(&data->A, c);
    CHOLMOD(finish)(c);

    qpdo_cleanup(work);
    c_free(data->q);
    c_free(data->l);
    c_free(data->u);
    c_free(data);
    c_free(settings);

    return 0;

}
