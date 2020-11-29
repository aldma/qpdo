# QPDO: the Quadratic Primal-Dual Optimizer

QPDO is a numerical solver for optimization problems in the form
```
minimize        0.5 x' Q x + q' x

subject to      l <= A x <= u
```
where `x in R^n` is the decision variable. The symmetric positive semidefinite matrix `Q in S_+^n`, the vector `q in R^n`, and the matrix `A in R^{m x n}` are bounded. The vectors `l in R^m U {-inf}^m` and `u in R^m U {+inf}^m` are extended-real-valued and satisfy `l_i ⩽ u_i` for all `i in 1,...,m`.


## Installation
QPDO is implemented in C and provides a MATLAB interface via mex.

Clone this repository with the submodule for SuiteSparse, running
```
git clone https://github.com/aldma/qpdo.git
git submodule update --init --recursive
```

### Matlab
* To install the mex interface of QPDO, add QPDO and its subfolders to the MATLAB path. Then go to [interfaces/mex/](./interfaces/mex/) and run `qpdo_make.m`. You can test and see how to call QPDO from MATLAB using `demo_mex.m` in the [examples/](./examples) folder.

## The Good and the Bad
Although this software package is still in its infancy, don't hesitate to [share with me](mailto:aldmarchi@gmail.com) your impression! Would you like to collaborate to build better software? I'm on board! Reporting mistakes is also very useful, and all types of issues are welcome, including bug reports, typos, and feature requests.

## Benchmarks
We are currently benchmarking QPDO against [OSQP](https://github.com/oxfordcontrol/osqp) and [QPALM](https://github.com/Benny44/QPALM), looking forward to sharing the results.
