function pb = init_problem(dim)
%INIT_PROBLEM  Specify the SPH problem
%   pb = init_problem(dim) initializes a problem for the specified
%   dimension D (1, 2, or 3).  The return value is a structure with the
%   following fields:
%      pb.dim     The problem dimension D
%      pb.domain  The ...
%      pb.n
%      pb.BC
%      ...

% Set spatial scale.  This is used for setting initial particle locations
% and setting the kernel's support radius.
h = 1;%%2e-4;

pb.dim = dim;

% density and pressure
pb.rho = 1;%%1180;
pb.pres = 1;%%1e8;

% initial particle velocity
% pb.v0 = zeros(pb.dim, 1);
% pb.v0(1) = 1;
pb.v0 = rand(pb.dim, 1);
pb.v0 = pb.v0 / norm(pb.v0);

% viscosity
pb.mu = 0.05;

% kernel support radius = 2 * pb.h
pb.h = 1.1 * h;
pb.eta2 = (0.1 * pb.h)^2;

% body force
pb.F = zeros(pb.dim, 1);
pb.F(1) = 1;%%4;

% time stepsize
pb.dt = 0.001;

% number of particles in each direction
% nx = 10;
% ny = 10;
% nz = 10;
nx = 5;
ny = 5;
nz = 5;

% domain limits
x_lim = [-3*h  3*h];
y_lim = [-3*h  3*h];
z_lim = [-3*h  3*h];

% types of BC in each direction.
% 1 - periodic
% 2 - wall
x_BC = 1;
y_BC = 2;
z_BC = 1;

switch pb.dim
    case 1
        pb.domain = x_lim;
        pb.BC = x_BC;
        pb.n = nx;
    case 2
        pb.domain = [x_lim; y_lim];
        pb.BC = [x_BC; y_BC];
        pb.n = [nx; ny];
        
    case 3
        pb.domain = [x_lim; y_lim; z_lim];
        pb.BC = [x_BC; y_BC; z_BC];
        pb.n = [nx; ny; nz];
end


% Particle mass
vol = prod(pb.domain(:,2) - pb.domain(:,1));
pb.m = pb.rho * vol / prod(pb.n);
