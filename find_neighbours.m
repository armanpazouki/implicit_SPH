function [nb_p, nb_g] = find_neighbours(loc, pb, part, ghost)
%FIND_NEIGHBOURS  Neighbour particles and ghosts to a given location
%   [nb_p, nb_g] = find_neighbours(loc, pb, part, ghost)
%   returns the index subsets nb_p (in 'part') and nb_g (in 'ghost') of
%   particles and ghosts, respectively, that are within the kernel support
%   radius of the specified point 'loc'.
nb_p = nearestneighbour(loc, part.r, 'r', 2.1 * pb.h);
nb_g = nearestneighbour(loc, ghost.r, 'r', 2.1 * pb.h);
