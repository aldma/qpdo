
clear all%#ok
close all
clc

% testing QPDO on degenerate and primal-dual infeasible problems

solver = qpdo;
settings                = solver.default_settings();
settings.verbose        = true;
settings.print_interval = 1;
settings.max_iter       = 100;

%% degenerate
qp_a = 0;
qp_b = 3;
qp_c = 0;

Q = [1,0; 0,0];
q = [1; qp_c];
A = [qp_a, qp_a; 1, 0; 0, 1];
l = [-inf; 1; 1];
u = [0; 3; qp_b];

solver = qpdo;
solver.setup(Q, q, A, l, u, settings);
res = solver.solve();
print_output( Q, q, A, l, u, res );

assert( res.info.status_val == 1 );

%% primal infeasible
qp_a = 1;
qp_b = 3;
qp_c = 0;

Q = [1,0; 0,0];
q = [1; qp_c];
A = [qp_a, qp_a; 1, 0; 0, 1];
l = [-inf; 1; 1];
u = [0; 3; qp_b];

solver = qpdo;
solver.setup(Q, q, A, l, u, settings);
res = solver.solve();
print_output( Q, q, A, l, u, res );

assert( res.info.status_val == -3 );

dy = res.prim_inf_cert;
norm_dy = norm( dy, inf );
disp( norm( A' * dy, inf ) / norm_dy );
out_of_bounds = dot( u(isfinite(u)), max(dy(isfinite(u)),0) ) + ...
                dot( l(isfinite(l)), min(dy(isfinite(l)),0) );
disp( out_of_bounds / norm_dy );

%% dual infeasible
% [0; 1] is an unbounded direction
% fprintf('QP is unbounded below, eflag should be -3 and dx should be a unbounded descent direction\n');
qp_a = 0;
qp_b = +inf;
qp_c = -1;

Q = [1,0; 0,0];
q = [1; qp_c];
A = [qp_a, qp_a; 1, 0; 0, 1];
l = [-inf; 1; 1];
u = [0; 3; qp_b];

solver = qpdo;
solver.setup(Q, q, A, l, u, settings);
res = solver.solve();
print_output( Q, q, A, l, u, res );

assert( res.info.status_val == -4 );

dx = res.dual_inf_cert;

norm_dx = norm( dx, inf );
disp( norm( Q * dx, inf ) / norm_dx );
disp( dot(q, dx) / norm_dx );
Adx = A*dx;
idx = isfinite(u) & isfinite(l);
disp( max(abs( Adx(idx)/norm_dx )) );

idx = isfinite(u) & isinf(l);
disp( max(Adx(idx)/norm_dx) );

idx = isinf(u) & isfinite(l);
disp( min(Adx(idx)/norm_dx) );



fprintf('That`s all folks! \n')

function print_output(Q, q, A, l, u, res)
    fprintf('%20s: %s\n', 'status',res.info.status);
    fprintf('%20s: %d\n', 'prox iterations',res.info.oterations);
    fprintf('%20s: %d\n', 'Newton iterations',res.info.iterations);
    if ~(any(isnan(res.x)) || any(isnan(res.y)))
        res_prim_norm = norm( A*res.x - max( l ,min( u , A*res.x + res.y )) , inf );
        res_dual_norm = norm( Q*res.x + q + A'*res.y                        , inf );
        fprintf('%20s: %7.4e\n','norm prim residual', res_prim_norm);
        fprintf('%20s: %7.4e\n','norm dual residual', res_dual_norm);
    end
    return
end