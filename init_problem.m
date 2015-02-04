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

% Pressure gradient and body force (in x direction)
% pb.K = 2e-9
pb.K = 0;%0.5;%2e-9; % 2 * mu^2 * Re / (rho*(Ly/2)^3)
pb.F = 0.0;

% Suggested number of SPH markers.
N = 100;

% Channel half-width and length
Ly = 1;
Lx = 1;%2*Ly;

% Calculate actual number of SPH markers and channel length.
pb.ny = floor(sqrt(Ly*N/Lx));
pb.del = Ly/(pb.ny + 1);
pb.nx = floor(N/pb.ny);
pb.N = pb.nx * pb.ny;
pb.Ly = Ly;
pb.Lx = pb.del * (pb.nx + 1);
pb.Ne = 0;
pb.NT = pb.Ne + pb.N;

% Estimate Reynolds number
u_max = (pb.K + pb.rho * pb.F) * (pb.Ly/2)^2 / (2 * pb.mu);
Re = u_max * (pb.Ly) / pb.nu;
fprintf('u_max = %g\n', u_max);
fprintf('Re = %g\n', Re);

% Dirichlet boundary
pb.uT = .001;
pb.vT = 0;
pb.uR = 0;
pb.vR = 0;
pb.uB = 0;
pb.vB = 0;
pb.uL = 0;
pb.vL = 0;

% Calculate particle mass.
pb.m = pb.rho * pb.del^2;

% Set the kernel support radius (2*h) such that we average N_avg neighbors.
N_avg = 30;
pb.h = 0.5 * pb.del * sqrt(N_avg / pi);

pb.eta2 = (0.1 * pb.h)^2;

% Time stepsize
% pb.dt = 10000;
pb.dt = 10;%10000; .01 is for normalized system

