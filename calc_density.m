function rho = calc_density(loc, pb, part, ghost, varargin)
%CALC_DENSITY  Calculate the density at a specified location.
%   rho = calc_density(loc, pb, part, ghost)
%      Force a evaluation of neighbour particles and ghosts.
%   rho = calc_density(loc, pb, part, ghost, nb_p, nb_g)
%      Use the sets 'nb_p' and 'n_g' of particle and ghost neighbours.
%
%   The density is calculated as:
%         rho = sum_{j in nb} {m_j * W(loc - r_j)}

if nargin == 4
    [nb_p, nb_g] = find_neighbours(loc, pb, part, ghost);
end

% Accumulate contributions from neighbour particles and ghosts.
rho = 0;

for i = 1:length(nb_p)
    d = loc - part.r(:,nb_p(i));
%     if abs(d(2)) < 1e-5
%         norm(d)
%     end
    w = kernel(d, pb.h, 0);
    rho = rho + pb.m * w;
end

for i = 1:length(nb_g)
    d = loc - ghost.r(:,nb_g(i));
    w = kernel(d, pb.h, 0);
    rho = rho + pb.m * w;
end
