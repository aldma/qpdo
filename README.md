# QPDO: the Quadratic Primal-Dual Optimizer

QPDO is a numerical solver for optimization problems in the form
```
minimize        0.5 x' Q x + q' x

subject to      l <= A x <= u
```
where `x in R^n` is the decision variable. The symmetric positive semidefinite matrix `Q in S_+^n`, the vector `q in R^n`, and the matrix `A in R^{m x n}` are bounded. The vectors `l in R^m U {-inf}^m` and `u in R^m U {+inf}^m` are extended-real-valued and satisfy `l_i ⩽ u_i` for all `i in 1,...,m`.

## Method and Citing
QPDO implements a primal-dual Newton proximal method for convex quadratic programming. The proposed method can handle degenerate problems, provides a mechanism for infeasibility detection, and can exploit warm starting, while requiring only convexity. In particular, all linear systems are solvable by construction, independently from the problem data, and an exact linesearch can be performed. Details can be found in the [research paper](https://doi.org/10.1007/s10589-021-00342-y) mentioned below, which serves as a user manual for advanced users. If you use QPDO in your work, we kindly ask that you cite the following reference.
```
@article{demarchi2022qpdo,
	author		= {De~Marchi, Alberto},
	title       	= {On a primal-dual {N}ewton proximal method for convex quadratic programs},
	journal     	= {Computational Optimization and Applications},
	year        	= {2022},
	volume      	= {81},
	number	    	= {2},
	pages 		= {369--395},
	doi         	= {10.1007/s10589-021-00342-y},
}
```

## Installation
QPDO is implemented in C and provides a MATLAB interface via mex, inspired by [OSQP](https://github.com/osqp/osqp) and [QPALM](https://github.com/Benny44/QPALM).

Clone this repository with the submodule for SuiteSparse, running
```
git clone https://github.com/aldma/qpdo.git
cd qpdo
git submodule update --init --recursive
```

### Matlab
* To install the mex interface of QPDO, add QPDO and its subfolders to the MATLAB path. Then go to [interfaces/mex/](./interfaces/mex/) and run `qpdo_make.m`. You can test and see how to call QPDO from MATLAB using `demo_mex.m` in the [examples/](./examples) folder.

## Get in touch
Don't hesitate to [share](mailto:aldmarchi@gmail.com) your impression! Would you like to collaborate to build better software? Here we are!
