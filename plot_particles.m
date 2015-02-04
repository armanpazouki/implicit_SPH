function [] = plot_particles(hf, part, ghost, pb)
%PLOT_PARTICLES  Plot the particles and ghosts.

figure(hf);


switch pb.dim
    case 1
        x_lim = pb.domain;
        x_min = x_lim(1) - 2*pb.h;
        x_max = x_lim(2) + 2*pb.h;
        
        y_lim = [0 0];
        y_min = y_lim(1) - 2*pb.h;
        y_max = y_lim(2) + 2*pb.h;
        
        plot([x_min x_max], [0 0], 'k-');
        plot([x_lim(1); x_lim(1)], [y_min y_max], 'k:');
        plot([x_lim(2); x_lim(2)], [y_min y_max], 'k:');
       
        plot(part.r, zeros(1,part.num), 'o');
        plot(ghost.r, zeros(1,ghost.num), 'co');
        
        quiver(part.r,zeros(1,part.num),part.v,zeros(1,part.num),0.3,'g')
        
        for i = 1 : part.num
            print_index_1D(part.r(i), i, pb.h);
        end
        
        for i = 1 : ghost.num
            print_index_1D(ghost.r(i), ghost.idx(i), pb.h);
        end
        
        axis equal
        xlim([x_min x_max]);
        ylim([y_min y_max]);
        
        set(gca, 'yticklabel', []);
    case 2
        x_lim = pb.domain(1,:);
        y_lim = pb.domain(2,:);
        
        x_min = x_lim(1) - 2*pb.h;
        x_max = x_lim(2) + 2*pb.h;
        
        y_min = y_lim(1) - 2*pb.h;
        y_max = y_lim(2) + 2*pb.h;
        
        rectangle('position', ...
                [x_lim(1) y_lim(1) x_lim(2)-x_lim(1) y_lim(2)-y_lim(1)], ...
                'edgecolor', 'none', ...
                'facecolor', [0.9 0.9 0.9]);
        colors = ['k', 'r'];
        styles = [':', '-'];
        plot([x_lim(1); x_lim(1)], [y_min y_max], 'linestyle', styles(pb.BC(1)), 'color', colors(pb.BC(1)));
        plot([x_lim(2); x_lim(2)], [y_min y_max], 'linestyle', styles(pb.BC(1)), 'color', colors(pb.BC(1)));
        plot([x_min; x_max], [y_lim(1) y_lim(1)], 'linestyle', styles(pb.BC(2)), 'color', colors(pb.BC(2)));
        plot([x_min; x_max], [y_lim(2) y_lim(2)], 'linestyle', styles(pb.BC(2)), 'color', colors(pb.BC(2)));

        plot(part.r(1,:), part.r(2,:), 'o');
        plot(ghost.r(1,:), ghost.r(2,:), 'co');
        
        quiver(part.r(1,:),part.r(2,:),part.v(1,:),part.v(2,:),0.3,'g')
        
        for i = 1 : part.num
            print_index_2D(part.r(:,i), i, pb.h);
        end
        
        for i = 1 : ghost.num
            print_index_2D(ghost.r(:,i), ghost.idx(i), pb.h);
        end
        
        axis equal
        xlim([x_min x_max]);
        ylim([y_min y_max]);
        
    case 3
        x_lim = pb.domain(1,:);
        y_lim = pb.domain(2,:);
        z_lim = pb.domain(3,:);
        
        x_min = x_lim(1) - 2*pb.h;
        x_max = x_lim(2) + 2*pb.h;
        
        y_min = y_lim(1) - 2*pb.h;
        y_max = y_lim(2) + 2*pb.h;
        
        z_min = z_lim(1) - 2*pb.h;
        z_max = z_lim(2) + 2*pb.h;

        corners = [x_lim(1) y_lim(1) z_lim(1)
                   x_lim(1) y_lim(2) z_lim(1)
                   x_lim(2) y_lim(2) z_lim(1)
                   x_lim(2) y_lim(1) z_lim(1)
                   x_lim(1) y_lim(1) z_lim(2)
                   x_lim(1) y_lim(2) z_lim(2)
                   x_lim(2) y_lim(2) z_lim(2)
                   x_lim(2) y_lim(1) z_lim(2)];
        faces = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
        patch('Vertices', corners, 'Faces', faces, 'FaceColor', [0.5 0.5 0.5], 'facealpha',0.2);
        
        scatter3(part.r(1,:), part.r(2,:), part.r(3,:), 10, 'k.');
        scatter3(ghost.r(1,:), ghost.r(2,:), ghost.r(3,:), 10, 'c.');
        
        view(40,35)
        axis equal
        xlim([x_min x_max]);
        ylim([y_min y_max]);
        zlim([z_min z_max]);
end

return

% =========================================================================

function [] = print_index_1D(loc, idx, del)
%PRINT_INDEX_1D  Print the specified index next to the given location.

loc = loc + del/15;
str = sprintf('%i', idx);
text(loc(1), del/15, str, ...
    'fontsize', 6, ...
    'horizontalalignment', 'left', ...
    'verticalalignment', 'bottom');

return

% =========================================================================

function [] = print_index_2D(loc, idx, del)
%PRINT_INDEX_2D  Print the specified index next to the given location.

loc = loc + del/15;
str = sprintf('%i', idx);
text(loc(1), loc(2), str, ...
    'fontsize', 6, ...
    'horizontalalignment', 'left', ...
    'verticalalignment', 'bottom');

return
