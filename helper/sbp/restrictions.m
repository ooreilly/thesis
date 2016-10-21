function r = restrictions(n)
% [e0p e0m eNp eNm] = restrictions(n)
% Construct the restrictions operators. For example,
% e0p = [1 0 ... 0]^T. 
% 
% Input arguments:
% n : number of grid points
% Output arguments
% r.e0p,r.eNp: Restriction operators for + grid (n+1) grid points
% r.e0m,r.eNm: Restriction operators for - grid (n+2) grid points

r.e0p = spalloc(n+1,1,1); r.e0p(1)   = 1;
r.eNp = spalloc(n+1,1,1); r.eNp(end) = 1;
r.e0m = spalloc(n+2,1,1); r.e0m(1)   = 1;
r.eNm = spalloc(n+2,1,1); r.eNm(end) = 1;

end
