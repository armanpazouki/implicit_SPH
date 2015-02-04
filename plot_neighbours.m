function plot_neighbours(hf, loc, pb, part, ghost)
%PLOT_NEIGHBOURS   Plot neighbours (particles and ghosts)

figure(hf);

% Find neighbour particles and ghosts for the specified location.
[nb_p, nb_g] = find_neighbours(loc, pb, part, ghost);

radius = 2 * pb.h;

switch pb.dim
    
    case 1
        h = rectangle('position',[loc-radius -pb.h/2 2*radius pb.h], 'curvature', [1 1]);
        set(h, 'linestyle', ':');
        
        plot(part.r(nb_p), zeros(size(nb_p)), 'o', 'markerfacecolor', 'y');
        plot(ghost.r(nb_g), zeros(size(nb_g)), 'co', 'markerfacecolor', 'y');
        
    case 2
        
        h = rectangle('position',[loc(1)-radius loc(2)-radius 2*radius 2*radius], 'curvature', [1 1]);
        set(h, 'linestyle', ':')
        
        plot(part.r(1,nb_p), part.r(2,nb_p), 'o', 'markerfacecolor', 'y');
        plot(ghost.r(1,nb_g), ghost.r(2,nb_g), 'co', 'markerfacecolor', 'y');
        
    case 3
        
        scatter3(part.r(1,nb_p), part.r(2,nb_p), part.r(3,nb_p), 30, 'o', 'markerfacecolor', 'y');
        scatter3(ghost.r(1,nb_g), ghost.r(2,nb_g), ghost.r(3,nb_g), 30, 'co', 'markerfacecolor', 'y');
        
        [x,y,z] = sphere(20);
        x = loc(1) + radius * x;
        y = loc(2) + radius * y;
        z = loc(3) + radius * z;
        
        h = surfl(x,y,z);
        set(h, 'FaceAlpha', 0.5, 'FaceColor', 'none','EdgeColor',0.85*[1 1 1]);
        
end
