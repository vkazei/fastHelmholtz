function D = F(m,model)
% Forward operator
%
%   D = P^TA^{-1}(m)Q
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   D = F(m,model);
%
% input:
%   m - squared-slownes [s^2/km^2]
%   model.h - gridspacing in each direction d = [d1, d2];
%   model.n - number of gridpoints in each direction n = [n1, n2]
%   model.f - frequency [Hz].
%   model.{zr,xr} - {z,x} locations of receivers [m] (must coincide with gridpoints)
%   model.{zs,xs} - {z,x} locations of sources [m] (must coincide with gridpoints)
%
%
% output:
%   D - data matrix
%
% Modified from https://github.com/TristanvanLeeuwen/SimpleFWI

% generate matrices
A = getA(model.f,m,model.h,model.n);
P = getP(model.zr,model.xr, model.z, model.x);
Q = getQ(model.zs,model.xs, model.z, model.x);

% solve
D = P'*(A\Q);
