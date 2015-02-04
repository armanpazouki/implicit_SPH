function part = init_particles(pb)
%INIT_PARTICLES  Initialize particle locations.
%    p = init_particles(pb) sets the initial particle locations based on
%    the problem specification 'pb'. The return value is a structure with
%    the following fields:
%       num    total number of particles, N
%       r      particle locations (2 x N matrix)
%       v      particle velocities (2 x N matrix)
%       p      particle pressures (1 x N vector)
%       nb_p   neighbour particles (N x 1 cell of arrays)
%       nb_g   neighbour ghosts (N x 1 cell of arrays)

%% Set initial particle locations
x = pb.del * (1:pb.nx);
y = pb.del * (1:pb.ny);
[X,Y] = meshgrid(x,y);
X = reshape(X, 1, pb.N);
Y = reshape(Y, 1, pb.N);
% part.r = Perturb([X;Y], 0.1*pb.del);
part.r = Perturb([X;Y], 0);


% Set initial particle velocities and pressures.
part.v = zeros(2, pb.N);
part.p = zeros(1, pb.N);


