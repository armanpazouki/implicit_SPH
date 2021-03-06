function pb = init_problem
%INIT_PROBLEM  Specify the SPH problem
%   pb = init_problem initializes the 2D problem.  


% Water density and kinematic and dynamic viscosity
% pb.rho = 1000;          % Kg/m^3
% pb.mu = 1e-3;           % Kg/(m^2*s^2)
% pb.nu = pb.mu / pb.rho; % m/s^2
pb.rho = 1000;          % Kg/m^3
pb.mu = 1;           % Kg/(m^2*s^2)
pb.nu = pb.mu / pb.rho; % m/s^2

% Suggested number of SPH markers.
N = 100;

% Channel half-width and length
b = .5;
L = 2*b;

pb.F = 1.6e-5;
% Estimate Reynolds number
% u_max = (pb.K + pb.rho * pb.F) * b^2 / (2 * pb.mu);
u_max = (pb.rho * pb.F) * b^2 / (2 * pb.mu) %2D Poiseuille

Re = u_max / 2 * (2 * b) / pb.nu
fprintf('u_max = %g\n', u_max);
fprintf('Re = %g\n', Re);

% Pressure gradient and body force (in x direction)
% pb.K = 2e-9
%pb.K = 0;%0.5;%2e-9; % 2 * mu^2 * Re / (rho*b^3)
pb.K = 0;%pb.rho * pb.F;%2 * pb.mu * u_max / b^2; %2D Poiseulle, delta pressure between the 2 ends

% Calculate actual number of SPH markers and channel length.
pb.ny = floor(sqrt(2*b*N/L));
pb.del = 2*b/pb.ny;
pb.nx = floor(N/pb.ny);
pb.N = pb.nx * pb.ny;
pb.b = b;
pb.L = pb.del * pb.nx;

% Calculate particle mass.
pb.m = pb.rho * pb.del^2;

% Set the kernel support radius (2*h) such that we average N_avg neighbors.
N_avg = 30;
pb.h = 0.5 * pb.del * sqrt(N_avg / pi);

pb.eta2 = (0.1 * pb.h)^2;

% Time stepsize
% pb.dt = 10000;
pb.dt = 1;%0.01;%10000; .01 is for normalized system

