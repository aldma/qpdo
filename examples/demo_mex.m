
addpath('../interfaces/mex/')

rng( 123456 );

%% problem
n           = 200;
m           = 100;
density     = 0.1;
rcondition  = 1e-3;
Q = sprandsym( n, density, rcondition, 2 );
q = randn( n, 1 );
A = sprandn( m, n, density, rcondition );
l = - rand( m, 1 );
u = + rand( m, 1 );
fprintf('Problem ready \n')
fprintf('\n\n');

%% solver
solver = qpdo;
settings                = solver.default_settings();
settings.max_iter       = 200;
settings.eps_abs        = 1e-6;
settings.verbose        = true;
settings.print_interval = 1;
solver.setup(Q, q, A, l, u, settings);
fprintf('Solver ready \n')
fprintf('\n\n');

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
fprintf('\n\n');

%% warm-start

x = x + 1e-3 * randn( n, 1 );
y = y + 1e-3 * randn( m, 1 );
solver.warm_start( x, y );
fprintf('Solver warm-started \n')
fprintf('\n\n');

res = solver.solve();

fprintf('Run time: %f s\n', res.info.run_time);
fprintf('Status: %s\n', res.info.status);

solver.delete();

fprintf('That`s all folks! \n')




