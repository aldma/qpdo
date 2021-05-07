
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
% fprintf('Degenerate but feasible QP, eflag should be 0 and dx, dv should be near zero\n');
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


% DY = res.prim_inf_cert;
% 
% einf = 1e-6 * norm(DY, inf);
% out_of_bounds = dot( u(isfinite(u)), max(DY(isfinite(u)), 0) ) + dot( l(isfinite(l)), min(DY(isfinite(l)), 0) );
% c1 = norm(A' * DY, inf) <= einf;
% c2 = out_of_bounds <= - einf;


%% primal infeasible
% fprintf('QP is  infeasible, eflag should be -2 and dv should be nonzero\n');
qp_a = 1;
qp_b = 3;
qp_c = 0;
Q = [1,0; 0,0];
q = [1; qp_c];
A = [qp_a,qp_a; 1,0; 0,1];
l = [-inf; 1; 1];
u = [0; 3; qp_b];

solver = qpdo;
solver.setup(Q, q, A, l, u, settings);
res = solver.solve();
print_output( Q, q, A, l, u, res );

%% dual infeasible
% [0; 1] is an unbounded direction
% fprintf('QP is unbounded below, eflag should be -3 and dx should be a unbounded descent direction\n');
qp_a = 0;
qp_b = +inf;
qp_c = -1;
Q = [1,0; 0,0];
q = [1; qp_c];
A = [qp_a,qp_a; 1,0; 0,1];
l = [-inf; 1; 1];
u = [0; 3; qp_b];

solver = qpdo;
solver.setup(Q, q, A, l, u, settings);
res = solver.solve();
print_output( Q, q, A, l, u, res );




function print_output(Q, q, A, l, u, res)
    x = res.x;
    y = res.y;
    res_prim_norm = norm( A*x - max( l ,min( u , A*x + y )) , inf );
    res_dual_norm = norm( Q*x + q + A'*y                    , inf );

    fprintf('%20s: %s\n', 'status',res.info.status);
    fprintf('%20s: %d\n', 'prox iterations',res.info.iterations);
    fprintf('%20s: %d\n', 'Newton iterations',res.info.oterations);
    fprintf('%20s: %7.4e\n','norm residual', max(res_dual_norm, res_prim_norm));
    return
end