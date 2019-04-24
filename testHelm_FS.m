% This is a simple test for Helmholtz solver with second order absorbing
% boundary conditions (ABC). We compare 1st and 2nd order ABCs with analytical solution
%
% (c) Vladimir Kazei, Oleg Ovcharenko and Dmitry Kabanov (2019)
% inspired by https://github.com/TristanvanLeeuwen/SimpleFWI

% Set colormap
set(groot,'DefaultFigureColormap',rdbuMap())

% Define model, [m]
dx = 10;
xmax = 16000;
zmax = 3000;

% Grid dimensions
n(2) = round(xmax/dx) + 1;
n(1) = round(zmax/dx) + 1;

% Homogeneous velocity model, [km/s]
v = 2 * ones(n); 

% Receivers array, [m]
xr = 0:2*dx:xmax;
zr = 0.5*zmax*ones(1,length(xr));

% Source location, [m]
xs = 0.5 * xmax;
zs = 1500;
% Source frequency, [Hz]
f = 5.0;

% Points per wavelength
lambda_min = min(v(:))/f;
ppw = lambda_min * 1000 / dx;
disp(['Points per wavelength ',num2str(ppw)])

% Grid, [m]
h  = dx * [1 1];
z  = [0:n(1)-1] * h(1);
x  = [0:n(2)-1] * h(2);
[zz,xx] = ndgrid(z,x);

% Squared slowness, [s/km]^2
m = 1./v(:).^2;

%% CONVENTIONAL, 1st order boundaries
% 1st order Helmholtz matrix
A1 = getA_1st_FS(f,m,h,n);
% Project wavefield to receiver locations
P = getP(zr,xr,z,x);
% Project wavefield to source locations
Q = getQ(zs,xs,z,x);
% 2nd order Helmholtz matrix
A2 = getA(f,m,h,n,true);

% Forward modeling, Conventional wavefield
tic;
U1 = A1\Q;
toc;
% Conventional wavefield solution for all sources [nz, nx]
U1_2D = reshape(U1,n);

%% ANALYTICAL
% Distance from source to each point in the model
r = @(zz,xx)((zz-zs).^2+(xx-xs).^2).^0.5;
% Angular frequency
omega = 1e-3*2*pi*f;
% Wavenumber
K = (omega/v(1));

% For 3D case
% G3D = @(zz,xx)exp(1i*K.*r(zz,xx))./r(zz,xx);

% Analytical wavefield solution. Green's funciton
G2D =  @(zz,xx)0.25i*besselh(0,2,conj(K).*r(zz,xx));

G2D = G2D(zz,xx);

 
%% symmetric source
zs = -zs;
r = @(zz,xx)((zz-zs).^2+(xx-xs).^2).^0.5;
G2D2 =  @(zz,xx)0.25i*besselh(0,2,conj(K).*r(zz,xx));

G2D2 = G2D2(zz,xx);

G_2D = conj(G2D - G2D2);



%G_2D = fillmissing(G_2D,'pchip');

%% ACCURATE, 2nd order boundaries
% Forward modeling, Accurate wavefield
%spparms('spumoni',2);
tic;
U2 = A2\Q;
toc;
% Accurate wavefield solution for all sources [nz, nx]
U2_2D = reshape(U2,n);

%% PLOT
close all;
% Wavefield (circles)
figure;
subplot 511;
imagesc(real(U1_2D)); 
axis equal tight; colorbar;
title('1st order boundaries');
caxis([-0.1 0.1]);

% Plot analytic solution
subplot 512;
imagesc(real(G_2D));
axis equal tight; colorbar;
title('Analytic wavefield');
caxis([-0.1 0.1]);

% Plot accurate wavefield
subplot 513;
imagesc(real(U2_2D));
axis equal tight; colorbar;
title('2nd order boundaries')
caxis([-0.1 0.1]);

% Plot phase difference Conventional<->Analytics
subplot 514;
imagesc(angle(U1_2D./G_2D));
caxis([-pi/4 pi/4]);
axis equal tight; colorbar;
title('Phase 1st - Analytics')

% Plot phase difference Accurate<->Analytics
subplot 515;
imagesc(angle(U2_2D./G_2D));
caxis([-pi/4 pi/4]);
axis equal tight; colorbar;
title('Phase 2nd-Analytics')

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
drawnow;

% %% TEST efficiency of banded solver
% model.xs = xr;
% model.zs = zr;
% model.xr = xr;
% model.zr = zr;
% model.f = f;
% model.h = h;
% model.n = n;
% 
% disp('Timing');
% 
% spparms('default');
% fprintf('\nStandard solver running...\n')
% tic;
% D = F(m, model);
% toc;
% 
% spparms('bandden',0);
% fprintf('\nBanded solver running...\n')
% tic;
% D = F(m, model);
% toc;
% 




