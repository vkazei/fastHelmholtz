function A = getA(f, m, h, n, fs)
% Create 2D Helmholtz matrix with 2nd order Clayton-Enquist absorbing boundary
%
% use:
%   A = getA(f,m,h,n,fs)
%
% input:
%   f - frequency [Hz]
%   m - squared-slownes [s^2/km^2]
%   h - gridspacing in each direction d = [d1, d2];
%   n - number of gridpoints in each direction n = [n1, n2]
%   fs - free surface flag, bool, true or false
%
% output:
%   A - sparse matrix
%
% Based on https://github.com/TristanvanLeeuwen/SimpleFWI
%
% Vladimir Kazei, Oleg Ovcharenko, 2019

%%
defval('fs',false);

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
% 2nd derivative in z direction
L11 = -kron(spy2,D1'*D1);

spy1 = speye(n(1));
spy1([1 end],[1 end])=0;
% 2nd derivative in x direction
L22 = -kron(D2'*D2,spy1);

% Laplacian inside the domain
L = L11 + L22;

% boundary
L1 = -kron(speye(n(2)),D1'*D1);
L2 = -kron(D2'*D2,speye(n(1)));

% 2nd order derivative along the boundary
L_BOUND = L2 + L1 - L11 - L22;

% this is a trick to have half Helmholtz in the second variable
% Matrix a is not equal to 1 at the boundary
a = ones(n);
a(:,[1 end]) = .5; a([1 end],:) = .5; a = a(:);

% boundary points of the domain
BND = find(a~=1);

% derivative with respect to the normal of the boundary
L_n = 0*L;

% Matrix a_corners is non-zero at the corners
a_corners = zeros(n);
a_corners(1,1) = 1;         % top left
a_corners(1,end) = 1;       % top right
a_corners(end,1) = 1;       % bottom left
a_corners(end,end) = 1;     % bottom right

% find linear indices of corner points
CORNERS = find(a_corners==1);
% adjusting 2nd derivative in the corners (distance is sqrt(2) shorter)
L_BOUND(CORNERS, :) = 2 * L_BOUND(CORNERS, :);

% Normal 1st derivatives at the corners
L_corner = 0*L;

% upper left corner
% compute the first derivative
L_corner(1,1) = 1;
L_corner(1,n(1)+2) = -1;
% distribute the second derivative
L_BOUND(1,1) = L_BOUND(1,1)/2;
L_BOUND(1,n(1)+2) = L_BOUND(1,1);

% lower left corner
% compute the first derivative
L_corner(n(1),n(1)) = 1;
L_corner(n(1),2*n(1) -1) = -1;
% distribute the second derivative
L_BOUND(n(1),n(1)) = L_BOUND(n(1),n(1))/2;
L_BOUND(n(1),2*n(1) -1) = L_BOUND(n(1),n(1));

% lower right corner
% compute the first derivative
L_corner(prod(n),prod(n)) = 1;
L_corner(prod(n),prod(n)-n(1)-1) = -1;
% distribute the second derivative
L_BOUND(prod(n),prod(n)) = L_BOUND(prod(n),prod(n))/2;
L_BOUND(prod(n),prod(n)-n(1)-1) = L_BOUND(prod(n),prod(n));


% upper right corner
% compute the first derivative
L_corner(prod(n)-n(1)+1, prod(n)-n(1)+1) = 1;
L_corner(prod(n)-n(1)+1, prod(n)-2*n(1)+2) = -1;
% distribute the second derivative
L_BOUND(prod(n)-n(1)+1, prod(n)-n(1)+1) =  L_BOUND(prod(n)-n(1)+1, prod(n)-n(1)+1)/2;
L_BOUND(prod(n)-n(1)+1, prod(n)-2*n(1)+2) =  L_BOUND(prod(n)-n(1)+1, prod(n)-n(1)+1);
%
L_corner = L_corner/(sqrt(2)*h(1));

% for the boundary points L is working as normal derivative
L_n(BND,:) = L(BND,:);

%% ASSEMBLE HELMHOLTZ MATRIX WITH 2ND ORDER ABC
% inside domain
A = diags(k.^2) + L-L_n;

% Clayton-Enquist 1977
% P_tt + (1/v) P_xt - (v/2) P_zz = 0
% we translate with the rules _t -> -i* omega;
% k^2 P + 1i * k P_x + 1/2 P_zz = 0
% where k = omega/v.

% Then we use prepared derivative operators
% _x -> L_n (normal derivative), _zz -> L_BOUND

A = A - 1i*(k(1))*(L_n * h(1) - L_corner) + 0.5*L_BOUND;

%A = prod(h)*A;

% Free surface
if fs
    b = zeros(n); b(1,:) = 1; b=b(:);
    FSP = find(b);
    A(FSP,:) = 0;
    A(:,FSP) = 0;
    A = omega^2*diags(b.*m) + A;
end

% check the number of grid points per wavelength
if (f > min(1e3*1./sqrt(m))/(5*h(1)))
    warning('Dispersion! f>min(1e3*1./sqrt(m))/(5*h(1))=%e',min(1e3*1./sqrt(m))/(5*h(1)));
end