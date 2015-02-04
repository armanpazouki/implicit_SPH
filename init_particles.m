function part = init_particles(pb)
%INIT_PARTICLES  Initialize particle locations.
%    p = init_particles(pb) sets the initial particle locations based on
%    the problem specification 'pb'. The return value is a structure with
%    the following fields:
%       num    total number of particles, N
%       r      particle locations (D x N matrix)
%       v      particle velocities (D x N matrix)
%       rho    particle densities (1 x N vector)
%       mu     particle viscosities (1 x N vector)
%       p      particle pressures (1 x N vector)
%       nb_p   neighbour particles (N x 1 cell of arrays)
%       nb_g   neighbour ghosts (N x 1 cell of arrays)

% Total number of particles.
part.num = prod(pb.n);

% Initialize particle densities, velocities, and pressures.
part.rho = pb.rho * ones(1, part.num);
part.mu  = pb.mu * ones(1, part.num);
part.v   = repmat(pb.v0, 1, part.num);
part.p   = pb.pres * ones(1, part.num);

% Initialize neighbour lists.
part.nb_p = cell(part.num, 1);
part.nb_g = cell(part.num, 1);

% Set initial particle locations
del = (pb.domain(:,2) - pb.domain(:,1))./ pb.n;

switch pb.dim
    
    case 1
        x = (pb.domain(1,1) + del(1)/2 : del(1) : pb.domain(1,2) - del(1)/2);
        part.r = x;
        
    case 2
        x = (pb.domain(1,1) + del(1)/2 : del(1) : pb.domain(1,2) - del(1)/2);
        y = (pb.domain(2,1) + del(2)/2 : del(2) : pb.domain(2,2) - del(2)/2);
        [X,Y] = meshgrid(x,y);
        X = reshape(X, 1, part.num);
        Y = reshape(Y, 1, part.num);
        part.r = [X;Y];
        
    case 3
        x = (pb.domain(1,1) + del(1)/2 : del(1) : pb.domain(1,2) - del(1)/2);
        y = (pb.domain(2,1) + del(2)/2 : del(2) : pb.domain(2,2) - del(2)/2);
        z = (pb.domain(3,1) + del(3)/2 : del(3) : pb.domain(3,2) - del(3)/2);
        [X,Y,Z] = meshgrid(x,y,z);
        X = reshape(X, 1, part.num);
        Y = reshape(Y, 1, part.num);
        Z = reshape(Z, 1, part.num);
        part.r = [X;Y;Z];
        
end

% m = -del/10;
% M = del/10;
% perturb = repmat(m, 1, part.num) + diag(M-m) * rand(size(part.r));
% part.r = part.r + perturb;

