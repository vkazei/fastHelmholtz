function A = getA_NEW(f,m,h,n)
% Create 2D Helmholtz matrix with 2nd order Clayton-Enquist absorbing boundary
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

%%
% scale frequency to handle velocity in km/s
omega = 1e-3*2*pi*f;

% z derivative operator
D1 = spdiags(ones(n(1),1)*[-1 1]/h(1),[0:1],n(1)-1,n(1));
% x derivative operator
D2 = spdiags(ones(n(2),1)*[-1 1]/h(2),[0:1],n(2)-1,n(2));

% internal points of the domain
spy2 = speye(n(2));
spy2([1 end],[1 end])=0;
L11 = -kron(spy2,D1'*D1);

spy1 = speye(n(1));
spy1([1 end],[1 end])=0;
L22 = -kron(D2'*D2,spy1);

L = L11+L22;

% boundary
L1 = -kron(speye(n(2)),D1'*D1);
L2 = -kron(D2'*D2,speye(n(1)));

L_BOUND = +0.5i*(h(1)^-1)*(L2+L1-L11-L22)/omega;

%this is a trick to have half Helmholtz in the second variable
a = ones(n); a(:,[1 end]) = .5; a([1 end],:) = .5; a = a(:);

tic;
BND = find(a~=1); 
mq = m.^0.25;
L_BOUND(BND,:) = L_BOUND(BND,:)./mq(BND);
L_BOUND(:,BND) = L_BOUND(:,BND)./((mq(BND))');
toc;

A = omega^2*diags(2*(a-0.5).*m) + (2i*omega/h(1))*diags((1-a).*sqrt(m)) + L + L_BOUND;

% check the number of grid points per wavelength
if (f > min(1e3*1./sqrt(m))/(5*h(1)))
    warning('Dispersion! f>min(1e3*1./sqrt(m))/(5*h(1))=%e',min(1e3*1./sqrt(m))/(5*h(1)));
end