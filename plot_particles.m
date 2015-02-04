function [] = plot_particles(hf, part, ghost, pb)
%PLOT_PARTICLES  Plot the particles and ghosts.

figure(hf);


x_lim = [0 pb.Lx];
y_lim = [-(0.5*pb.Ly) (0.5*pb.Ly)];

x_min = x_lim(1) - 2*pb.h;
x_max = x_lim(2) + 2*pb.h;

y_min = y_lim(1) - 2*pb.h;
y_max = y_lim(2) + 2*pb.h;

rectangle('position', ...
    [x_lim(1) y_lim(1) x_lim(2)-x_lim(1) y_lim(2)-y_lim(1)], ...
    'edgecolor', 'none', ...
    'facecolor', [0.9 0.9 0.9]);
plot([x_lim(1); x_lim(1)], [y_min y_max], 'k-');
plot([x_lim(2); x_lim(2)], [y_min y_max], 'k-');
plot([x_min; x_max], [y_lim(1) y_lim(1)], 'r:');
plot([x_min; x_max], [y_lim(2) y_lim(2)], 'r:');

plot(part.r(1,:), part.r(2,:), 'o');
plot(ghost.r(1,:), ghost.r(2,:), 'co');

quiver(part.r(1,:),part.r(2,:),part.v(1,:),part.v(2,:),0.3,'g')

for i = 1 : pb.N
    print_index_2D(part.r(:,i), i, pb.del);
end

for i = 1 : ghost.num
    print_index_2D(ghost.r(:,i), ghost.idx(i), pb.del);
end

axis equal
xlim([x_min x_max]);
ylim([y_min y_max]);

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
