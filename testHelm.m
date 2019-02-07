%% this is a simple test for Helmholtz
set(groot,'DefaultFigureColormap',rdbuMap())

% Define model
dx = 40;
x_length = 16000;
z_length = 16000;

% Grid dimensions
n(2) = round(x_length/dx)+1;
n(1) = round(z_length/dx)+1;

% Homogeneous velocity model 2 km/s
v = 2 * ones(n);

% Receivers array
xr = 0:dx:2000;
zr = dx*ones(1,length(xr));

% Source
xs = 2000;      % location
zs = 2000;
f = 1.0;        % source frequncy

% Points per wavelength
lambda_min = min(v(:))/f;
ppw = lambda_min * 1000 / dx;

% Grid
h  = dx * [1 1];
z  = [0:n(1)-1] * h(1);
x  = [0:n(2)-1] * h(2);
[zz,xx] = ndgrid(z,x);

% Squared slowness
m = 1./v(:).^2;

%% 
% 1st order Helmholtz matrix
A1 = getA_1st(f,m,h,n);
% Project wavefield to receiver locations
P = getP(h,n,zr,xr);
% Project wavefield to source locations
Q = getP(h,n,zs,xs);
% 2nd order Helmholtz matrix
A2 = getA(f,m,h,n);

%% CONVENTIONAL, 1sr order boundaries
% Resource monitor
% Comment out the following line to mute the output
spparms('spumoni',2);

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
G_2D = @(zz,xx)0.25i * besselh(0,2,conj(K) .* r(zz,xx));
G_2D = conj(G_2D(zz,xx));

%% ACCURATE, 2nd order boundaries
% Forward modeling, Accurate wavefield
spparms('spumoni',2);
tic;
U2 = A2\Q;
toc;
U2_2D = reshape(U2,n);

%% PLOT
close all;
% Wavefield (circles)
figure;
subplot 321;
imagesc(real(U1_2D)); 
axis equal tight; colorbar;
title('1st order boundaries');
caxis([-0.1 0.1]);

% Plot analytic solution
subplot 322;
imagesc(real(G_2D));
axis equal tight; colorbar;
title('Analytic wavefield');
caxis([-0.1 0.1]);

% Plot phase difference Conventional<->Analytics
subplot 323;
imagesc(angle(U1_2D./G_2D));
caxis([-pi/4 pi/4]);
axis equal tight; colorbar;
title('Phase 1st - Analytics')

% Plot phase difference Accurate<->Analytics
subplot 325;
imagesc(angle(U2_2D./G_2D));
caxis([-pi/4 pi/4]);
axis equal tight; colorbar;
title('Phase 2nd-Analytics')

% Plot accurate wavefield
subplot 324;
imagesc(real(U2_2D));
axis equal tight; colorbar;
title('2nd order boundaries')
caxis([-0.1 0.1]);

% Plot relative amplitude error Accurate<->Analytics
subplot 326;
imagesc(1-abs(U2_2D./G_2D));
axis equal tight; colorbar;
title('Relative amplitude error 2nd-Analytics')
caxis([-0.2 0.2]);


% Plot slices from the middle of the model
figure;
% Conventional wavefield
plot(real(U1_2D(round(n(1)/2),:)),'g','linewidth',2); hold on;
% Accurate wavefield
plot(real(U2_2D(round(n(1)/2),:)),'b','linewidth',2);
% Analytic wavefield
plot(real(G_2D(round(n(1)/2),:)),'r','linewidth',2);
legend('1st order boundaries', '2nd order boundaries', 'Analytic');

