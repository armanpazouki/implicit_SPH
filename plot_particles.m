function [] = plot_particles(hf, part, edgeP, dummyP, pb)
%PLOT_PARTICLES  Plot the particles and ghosts.

figure(hf);


x_lim = [0 pb.Lx];
y_lim = [0 pb.Ly];

x_min = x_lim(1) - 4*pb.h;
x_max = x_lim(2) + 4*pb.h;

y_min = y_lim(1) - 4*pb.h;
y_max = y_lim(2) + 4*pb.h;

rectangle('position', ...
    [x_lim(1) y_lim(1) x_lim(2)-x_lim(1) y_lim(2)-y_lim(1)], ...
    'edgecolor', 'none', ...
    'facecolor', [0.9 0.9 0.9]);
plot([x_lim(1); x_lim(1)], [y_min y_max], 'k-');
plot([x_lim(2); x_lim(2)], [y_min y_max], 'k-');
plot([x_min; x_max], [y_lim(1) y_lim(1)], 'r:');
plot([x_min; x_max], [y_lim(2) y_lim(2)], 'r:');

plot(part.r(1,:), part.r(2,:), 'bo');
plot(edgeP.r(1,:), edgeP.r(2,:), 'co');
plot(dummyP.r(1,:), dummyP.r(2,:), 'ro');

quiver(part.r(1,:),part.r(2,:),part.v(1,:),part.v(2,:),0.3,'b')
quiver(edgeP.r(1,:),edgeP.r(2,:),edgeP.v(1,:),edgeP.v(2,:),0.3,'c')
quiver(dummyP.r(1,:),dummyP.r(2,:),edgeP.v(1, dummyP.idx),edgeP.v(2, dummyP.idx),0.3,'r')

for i = 1 : pb.N
    print_index_2D(part.r(:,i), i, pb.del);
end

for i = 1 : edgeP.num
    print_index_2D(edgeP.r(:,i), edgeP.idx(i), pb.del);
end

for i = 1 : dummyP.num
    print_index_2D(dummyP.r(:,i), dummyP.idx(i), pb.del);
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
