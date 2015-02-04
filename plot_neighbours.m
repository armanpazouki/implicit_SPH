function plot_neighbours(hf, loc, pb, part, ghost)
%PLOT_NEIGHBOURS   Plot neighbours (particles and ghosts)

figure(hf);

% Find neighbour particles and ghosts for the specified location.
[nb_p, nb_g] = find_neighbours(loc, pb, part.r, ghost.r);

radius = 2 * pb.h;

h = rectangle('position',[loc(1)-radius loc(2)-radius 2*radius 2*radius], 'curvature', [1 1]);
set(h, 'linestyle', ':')

plot(part.r(1,nb_p), part.r(2,nb_p), 'o', 'markerfacecolor', 'y');
plot(ghost.r(1,nb_g), ghost.r(2,nb_g), 'co', 'markerfacecolor', 'y');
        
