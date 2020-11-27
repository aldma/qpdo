# QPDO: the Quadratic Primal-Dual Optimizer

QPDO is a numerical solver for optimization problems in the form
```
minimize        0.5 x' Q x + q' x

subject to      l <= A x <= u
```
where `x in R^n` is the decision variable. The symmetric positive semidefinite matrix `Q in S_+^n`, the vector `q in R^n`, and the matrix `A in R^{m x n}` are bounded. The vectors `l in R^m U {-inf}^m` and `u in R^m U {+inf}^m` are extended-real-valued and satisfy `l_i ⩽ u_i` for all `i in 1,...,m`.
