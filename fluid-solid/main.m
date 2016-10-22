clear all
addpath('../helper/sbp');
addpath('../helper/sparse');
addpath('../helper/tests');
addpath('../helper/time_integrators');
addpath('../helper/solvers');
addpath('../helper/pgfplots');
addpath('../pde/');

n     = 300;
order = 4;
CFL   = 0.5;
n1    = n;
n2    = n;
tend  = 1;
show_plot = true;
frame_stride = 5;
save_at_time = 5;

f = Fluid(n,order);

f.test_stability(false)

s = Solver(n,n,order);
s = s.solve(CFL,tend,show_plot,frame_stride,save_at_time);
