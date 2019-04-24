function P = getP(zs,xs,z,x)
% Define sampling operator
%
% use:
%   P = getP(h,n,zt,xt)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%   zs,xs - arrays defining sampling points (must coincide with grid)
%
% output
%   P     - sparse matrix
%
% Modified from https://github.com/TristanvanLeeuwen/SimpleFWI


[zz,xx] = ndgrid(z,x);

for k = 1:length(zs)
    i(k) = find((zz(:)==zs(k))&(xx(:)==xs(k)));
end

I = speye(length(x)*length(z));

P = I(:,i);
