function [nb_p, nb_d] = find_neighbours(loc, pb, r_part_total, r_dummyP)
%FIND_NEIGHBOURS  Neighbour particles and ghosts to a given location
%   [nb_p, nb_g] = find_neighbours(loc, pb, part, ghost)
%   returns the index subsets nb_p (in 'part') and nb_g (in 'ghost') of
%   particles and ghosts, respectively, that are within the kernel support
%   radius of the specified point 'loc'.
nb_p = nearestneighbour(loc, r_part_total, 'r', 2.1 * pb.h);
nb_d = nearestneighbour(loc, r_dummyP, 'r', 2.1 * pb.h);
