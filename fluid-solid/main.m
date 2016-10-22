clear all
addpath('../helper/sbp');
addpath('../helper/sparse');
addpath('../helper/tests');
addpath('../helper/time_integrators');
addpath('../helper/solvers');
addpath('../helper/pgfplots');
addpath('../pde/');

n     = 300;
order = 6;
CFL   = 0.5;
n1    = 100*n;
n2    = n;
tend  = 1;
show_plot = true;
frame_stride = 5;
save_at_time = [0 0.3 0.6 0.9];

fluid = Fluid(n,order);
fluid.mu = 1e-3;
solid = Solid(n,order);


s = Solver(fluid,solid);
s = s.solve(CFL,tend,show_plot,frame_stride,save_at_time);
