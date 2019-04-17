function A = getA_1st_FS(f,m,h,n)
% Define 2D Helmholtz matrix with absorbing boundary conditions
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
% Upgraded to free-surface case from https://github.com/TristanvanLeeuwen/SimpleFWI
%
% Vladimir Kazei, Oleg Ovcharenko, 2019

%%
omega = 1e-3*2*pi*f;
N     = prod(n);

D1 = spdiags(ones(n(1),1)*[-1 1]/h(1),[0:1],n(1)-1,n(1));
D2 = spdiags(ones(n(2),1)*[-1 1]/h(2),[0:1],n(2)-1,n(2));
spy2 = speye(n(2));
spy2([1 end],[1 end])=0;
L11 = -kron(spy2,D1'*D1);
L1 = -kron(speye(n(2)),D1'*D1);
spy1 = speye(n(1));
spy1([1 end],[1 end])=0;
L22 = -kron(D2'*D2,spy1);
L2 = -kron(D2'*D2,speye(n(1)));
L = L11+L22;

L_BOUND = +0.5i*(h(1)^-1)*(L2+L1-L11-L22)/omega;

%this is a trick to have half Helmholtz in the second variable
a = ones(n); a(:,[1 end]) = .5; a([1 end],:) = .5; a = a(:);

BND = find(a~=1); 
mq = m.^0.5;

L_BOUND(BND,:) = L_BOUND(BND,:)./mq(BND);
L_BOUND(:,BND) = L_BOUND(:,BND)./((mq(BND))');

A = (2i*omega/h(1))*diags((1-a).*sqrt(m)) + L + L_BOUND;

%% free surface
b = zeros(n); b(1,:) = 1; b=b(:);
FSP = find(b);

% zero out the elements out of the diagonals in A for free surface elements
A(FSP,:) = 0;
A(:,FSP) = 0;
A = omega^2*diags(2*(a-0.5).*m) + omega^2*diags(b.*m) + A;
%A = prod(h)*A;

if (f > min(1e3*1./sqrt(m))/(7.5*h(1)))
    warning('f>min(1e3*1./sqrt(m))/(7.5*h(1))=%e',min(1e3*1./sqrt(m))/(7.5*h(1)));
end