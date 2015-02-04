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

%% Initialize with fully-developped flow?
fully_dev = false;%true;


% Initialize neighbour lists.
part.nb_p = cell(pb.N, 1);
part.nb_g = cell(pb.N, 1);

% Set initial particle locations
x = pb.del * (1:pb.nx) - pb.del/2;
y = pb.del * (1:pb.ny) - pb.b - pb.del/2;
[X,Y] = meshgrid(x,y);
X = reshape(X, 1, pb.N);
Y = reshape(Y, 1, pb.N);
% part.r = Perturb([X;Y], 0.1*pb.del);
part.r = Perturb([X;Y], 0);


% Set initial particle velocities and pressures.
part.v = zeros(2, pb.N);
part.p = zeros(1, pb.N);

if fully_dev
    K = pb.K + pb.rho * pb.F;
    tmp = K * pb.b^2 / (2 * pb.mu);
    
    for i = 1:pb.N
        part.v(1,i) = tmp * (1-(part.r(2,i)/pb.b)^2);
%         part.p(i) = pb.K * (pb.L - part.r(1,i));
    end
    
end


% m = -pb.del/10;
% M = pb.del/10;
% perturb = repmat(m, 1, part.num) + diag(M-m) * rand(size(part.r));
% part.r = part.r + perturb;

