function [part, ghost] = pre_process(part, pb)
%

%PRE_PROCESS This function performs the step pre-processing:
%   - set the locations and associated particle for the ghosts.
%   - calculate the set of neighbours (other particles and/or ghosts) for
%     each particle.
%   - calculate the density at each particle location.
%   - set the density and pressure at each ghost location.

ghost = set_ghosts(part, pb);

for i = 1 : part.num
    [nb_p, nb_g] = find_neighbours(part.r(:,i), pb, part, ghost);
    part.nb_p{i} = nb_p;
    part.nb_g{i} = nb_g;
    part.rho(i) = calc_density(part.r(:,i), part, ghost, pb, nb_p, nb_g);
end


ghost.rho = part.rho(ghost.idx);
ghost.p = part.p(ghost.idx);
ghost.v = part.v(ghost.idx);
