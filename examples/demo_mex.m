
addpath('../interfaces/mex/')

%% problem
n           = 30;
m           = 50;
density     = 0.1;
rcondition  = 1e-2;
Q = sprandsym( n, density, rcondition, 2 );
q = randn( n, 1 );
A = sprandn( m, n, density, rcondition );
l = - rand( m, 1 );
u = + rand( m, 1 );
fprintf('Problem ready!\n')

%% solver
solver = qpdo;
settings                = solver.default_settings();
settings.max_iter       = 200;
settings.eps_abs        = 1e-6;
settings.verbose        = true;
settings.print_interval = 1;
solver.setup(Q, q, A, l, u, settings);
fprintf('Solver ready!\n')

%% cold-start
res = solver.solve();

fprintf('Run time: %f s\n', res.info.run_time);
fprintf('Status: %s\n', res.info.status);

x = res.x;
y = res.y;

res_prim_norm = norm( A*x - max( l ,min( u , A*x + y )) , inf );
res_dual_norm = norm( Q*x + q + A'*y                    , inf );
fprintf('primal res: %7.4e\n', res_prim_norm);
fprintf('  dual res: %7.4e\n', res_dual_norm);

%% warm-start
x = x + randn( n, 1 );
y = y + randn( m, 1 );
solver.warm_start( x, y );
fprintf('Solver ready!\n')

res = solver.solve();

fprintf('Run time: %f s\n', res.info.run_time);
fprintf('Status: %s\n', res.info.status);

solver.delete();
fprintf('Solver deleted!\n')




