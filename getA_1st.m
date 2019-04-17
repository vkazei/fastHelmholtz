function A = getA_1st(f,m,h,n)
% Create 2D Helmholtz matrix with 1st order Enquist absorbing boundary
%
% use:
%   A = getA(f,m,h,n)
%
% input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   h - gridspacing in each direction d = [d1, d2];
%   n - number of gridpoints in each direction n = [n1, n2]
%
% output:
%   A - sparse matrix
%
% Modified from https://github.com/TristanvanLeeuwen/SimpleFWI

%%
% scale frequency to handle velocity in km/s
omega = 1e-3*2*pi*f;
k = omega.*sqrt(m);

% z derivative operator
D1 = spdiags(ones(n(1),1)*[1 -1]/h(1),[0:1],n(1)-1,n(1));
% x derivative operator
D2 = spdiags(ones(n(2),1)*[1 -1]/h(2),[0:1],n(2)-1,n(2));

% internal points of the domain
spy2 = speye(n(2));
% zeroing out the boundary
spy2([1 end],[1 end])=0;
L11 = -kron(spy2,D1'*D1);

spy1 = speye(n(1));
spy1([1 end],[1 end])=0;
L22 = -kron(D2'*D2,spy1);

% Laplacian inside 
L = L11+L22;

% this is a trick to have half Helmholtz in the second variable
a = ones(n); a(:,[1 end]) = .5; a([1 end],:) = .5; a = a(:);

% boundary points of the domain
BND = find(a~=1);

% normal derivative operator
L_n = 0*L;

L_corner = 0*L;

% upper left corner
L_corner(1,1) = 1;
L_corner(1,n(1)+1) = -1;

% lower left corner 
L_corner(n(1),n(1)) = 1;
L_corner(n(1),2*n(1) -1) = -1;

% lower right corner
L_corner(prod(n),prod(n)) = 1;
L_corner(prod(n),prod(n)-n(1)-1) = -1;

% upper right corner
L_corner(prod(n)-n(1)+1, prod(n)-n(1)+1) = 1;
L_corner(prod(n)-n(1)+1, prod(n)-2*n(1)+2) = -1;

L_corner = L_corner/(sqrt(2)*h(1));

% for the boundary points L is working as normal derivative
L_n(BND,:) = L(BND,:);

% inside domain
A = diags(k.^2) + L-L_n;

% Clayton-Enquist 1977
% P_t + (1/v) P_x = 0
% we translate with the rules _t -> -i* omega; 
% k^2 P + 1i * k P_x = 0 
% where k = omega/v.

% Then we use prepared derivative operators
% _x -> L_n (normal derivative), _zz -> L_BOUND

A = A - 1i*(k(1))*(L_n * h(1) - L_corner);

%A = prod(h)*A;

% check the number of grid points per wavelength
if (f > min(1e3*1./sqrt(m))/(5*h(1)))
    warning('Dispersion! f>min(1e3*1./sqrt(m))/(5*h(1))=%e',min(1e3*1./sqrt(m))/(5*h(1)));
end