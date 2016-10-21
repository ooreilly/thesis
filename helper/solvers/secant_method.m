function xs  = secant_method(fun,x0,xl,xr,tol_x,tol_f,debug)
% secant_method(fun,xl,xr,tol_x,debug)
% Finds a root to fun(x) = 0.
%
% Input:
%    fun: fun = @(dt) ... Function handle to stability function to minimize
%    x0 : initial guess for center
%    xl : initial guess for left bound
%    xr : initial guess for right bound
% Optional:
%    tol_x: Tolerance for the stopping criterion  xr - xl < tol.
%    debug: Display iteration information
% Output:
%     xs: The root

if nargin < 4
  tol_x = 1e-3;
end
if nargin < 5
  tol_x = 1e-6;
end
if nargin < 6
  debug = false;
end

f0 = fun(x0);
fl = fun(xl);
fr = fun(xr);

if f0 == 0
  xs = x0;
  return;
end

if fl == 0
  xs = xl;
  return;
end
if fr == 0
  xs = xr;
  return;
end

dx = xr - xl;

df = abs(f0);
i = 1;
  
x0 = (xl+xr)/2;
if debug
  fprintf('Iter. \t  lb \t ub \t err f \n');
end

while dx > tol_x || df > tol_f 

  fl = fun(xl);
  fr = fun(xr);
  f0 = fun(x0);

  search_left  = root_in_bounds(fl,f0);
  search_right = root_in_bounds(f0,fr);

  if search_left
    xl = xl;
    xr = x0;
  end
  if search_right
    xl = x0;
    xr = xr;
  end

  % Update bounds
  dx = xr - xl;

  df = abs(f0);

  if debug
    fprintf('%05d \t %.2g \t %.2g \t %.2g \n',i,xl,xr,abs(f0));
  end
  i = i + 1;
  
  % Split interval into two intervals: (xl,x0) (x0,xr)
  x0 = (xl+xr)/2;

end

xs = x0;


end

function y = root_in_bounds(fl,fr)
  tol_f = 1e-3;
  if fl*fr < 0
    y = true;
  else
    y = false;
  end
end
  
