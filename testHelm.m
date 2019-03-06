% This is a simple test for Helmholtz solver with second order absorbing
% boundary conditions (ABC). We compare 1st and 2nd order ABCs with analytical solution
% (c) Vladimir Kazei, Oleg Ovcharenko and Dmitry Kabanov (2019)
% inspired by https://github.com/TristanvanLeeuwen/SimpleFWI

set(groot,'DefaultFigureColormap',rdbuMap())

% Define model
dx = 20;
x_length = 16000;
z_length = 4000;

% Grid dimensions
n(2) = round(x_length/dx)+1;
n(1) = round(z_length/dx)+1;

% Homogeneous velocity model 2 km/s
v = 2 * ones(n); 

% Receivers array
xr = 0:2*dx:x_length;
zr = 0.5*z_length*ones(1,length(xr));

% Source
xs = 0.25*x_length;     % location
zs = 2000;
f = 5;                % source frequency

% Points per wavelength
lambda_min = min(v(:))/f;
ppw = lambda_min * 1000 / dx;
disp(['Points per wavelength ',num2str(ppw)])

% Grid
h  = dx * [1 1];
z  = [0:n(1)-1] * h(1);
x  = [0:n(2)-1] * h(2);
[zz,xx] = ndgrid(z,x);

% Squared slowness
m = 1./v(:).^2;

%% CONVENTIONAL, 1st order boundaries
% 1st order Helmholtz matrix
A1 = getA_1st(f,m,h,n);
% Project wavefield to receiver locations
P = getP(h,n,zr,xr);
% Project wavefield to source locations
Q = getP(h,n,zs,xs);
% 2nd order Helmholtz matrix
A2 = getA(f,m,h,n);

% Forward modeling, Conventional wavefield
tic;
U1 = A1\Q;
toc;

% Conventional wavefield solution for all sources [nz, nx]
U1_2D = reshape(U1,n);

%% ANALYTICAL
% Distance from source to each point in the model
r = @(zz,xx)(zz.^2+xx.^2).^0.5;
% Angular frequency
omega = 1e-3*2*pi*f;
% Wavenumber
K = (omega/v(1));

% For 3D case
% G3D = @(zz,xx)exp(1i*K.*r(zz,xx))./r(zz,xx);

% Analytical wavefield solution. Green's funciton
G_2D_analytic = @(zz,xx)0.25i * besselh(0,2,conj(K) .* r(zz,xx));
G_2D = conj(G_2D_analytic(zz - zs, xx - xs));

%G_2D = fillmissing(G_2D,'pchip');

%% ACCURATE, 2nd order boundaries
% Forward modeling, Accurate wavefield
%spparms('spumoni',2);
tic;
U2 = A2\Q;
toc;
U2_2D = reshape(U2,n);

%% PLOT
close all;
% Wavefield (circles)
figure;
subplot 231;
imagesc(real(U1_2D)); 
axis equal tight; colorbar;
title('1st order boundaries');
caxis([-0.1 0.1]);

% Plot analytic solution
subplot 232;
imagesc(real(G_2D));
axis equal tight; colorbar;
title('Analytic wavefield');
caxis([-0.1 0.1]);

% Plot accurate wavefield
subplot 233;
imagesc(real(U2_2D));
axis equal tight; colorbar;
title('2nd order boundaries')
caxis([-0.1 0.1]);

% Plot phase difference Conventional<->Analytics
subplot 234;
imagesc(angle(U1_2D./G_2D));
caxis([-pi/4 pi/4]);
axis equal tight; colorbar;
title('Phase 1st - Analytics')

% Plot phase difference Accurate<->Analytics
subplot 235;
imagesc(angle(U2_2D./G_2D));
caxis([-pi/4 pi/4]);
axis equal tight; colorbar;
title('Phase 2nd-Analytics')

% Plot relative amplitude error Accurate<->Analytics
subplot 236;
imagesc(1-abs(U2_2D./G_2D));
axis equal tight; colorbar;
title('Relative amplitude error 2nd-Analytics')
caxis([-0.2 0.2]);

% Print out relative improvement, %
diff1 = fillmissing(G_2D-U1_2D, 'linear');
diff2 = fillmissing(G_2D-U2_2D, 'linear');
disp('relative improvement of 2nd order absorbing boundaries over 1st, %')
disp(100*(norm(diff1, 'fro') - norm(diff2, 'fro')) / norm(diff1,'fro'))

%% CHECK HELMHOLTZ MATRIX ACTION ON TRUE GREEN'S FUNCITON
% Ideally, a delta function in source location should show up

% 1st order BC 
source_1_G2D = abs(A1*G_2D(:));
% 2nd order BC
source_2_G2D = abs(A2*G_2D(:));

% Plot 1st order product
figure;
subplot 211
title 1st
imagesc(reshape(source_1_G2D, size(G_2D)));
caxis([-1 1]*1e-3);
cax = caxis();

% Plot 2nd order product
subplot 212
title 2nd
imagesc(reshape(source_2_G2D, size(G_2D)));
caxis(cax);

%% TEST efficiency of banded solver
model.xs = xr;
model.zs = zr;
model.xr = xr;
model.zr = zr;
model.f = f;
model.h = h;
model.n = n;
spparms('default');
disp('Standard MATLAB solver')
tic;
D = F(m, model);
toc;
spparms('bandden',0);
disp('Banded solver')
tic;
D = F(m, model);
toc;





