clear
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

show_plot    = true;
frame_stride = 5;
tend         = 1.0;
save_at_times = [0.05 0.2 0.4 0.6 0.8];

f      = Fault();
f.taul = @(t)0.5*(1 - exp(-5*t));


s = Solver(CFL,n1,n2,order,f,tend);
s = s.solve(show_plot,frame_stride,save_at_times);
