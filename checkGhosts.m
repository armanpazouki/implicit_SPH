%CHECKGHOSTS  Check and plot the particles and ghosts for a 2-D problem.
%   The check consists of calculating the SPH density estimate at various
%   points and comparing it against the base value.

pb = init_problem();
part = init_particles(pb);

%% Set the ghost points.
ghost = set_ghosts(pb, part);

% Plot particles and ghosts
hf = figure;
set(hf, 'position', [50, 50, 640, 640]);
hold on

plot_particles(hf, part, ghost, pb);

%% Find and plot neighbors for a few locations.
% Estimate density at those locations.
fprintf('Specified density:  %f\n', pb.rho);

% center location.
loc = [0;0];
rho = calc_density(loc, pb, part, ghost);
fprintf('  at [%f %f]  rho = %f\n', loc(1), loc(2), rho);

% location close to right-bottom corner.
loc = [2*pb.h; -2.5*pb.h];
plot_neighbours(hf, loc, pb, part, ghost);
plot(loc(1), loc(2), 'r*', 'markerfacecolor', 'r');
rho = calc_density(loc, pb, part, ghost);
fprintf('  at [%f %f]  rho = %f\n', loc(1), loc(2), rho);

% particle close to top boundary.
ip = round(pb.nx/2);
jp = pb.ny-1;
loc = part.r(:,(ip-1)*pb.ny + jp);
plot_neighbours(hf, loc, pb, part, ghost);
plot(loc(1), loc(2), 'r*', 'markerfacecolor', 'r');
rho = calc_density(loc, pb, part, ghost);
fprintf('  at [%f %f]  rho = %f\n', loc(1), loc(2), rho);

%% Plot pressure
x_lim = [0 pb.Lx];
y_lim = [-(0.5*pb.Ly) (0.5*pb.Ly)];

x_min = x_lim(1) - 2*pb.h;
x_max = x_lim(2) + 2*pb.h;

y_min = y_lim(1) - 2*pb.h;
y_max = y_lim(2) + 2*pb.h;


hf = figure;
set(hf, 'position', [200, 80, 800, 600]);

subplot(3,1,1)
x = part.r(1,:)';
y = part.r(2,:)';
z = part.p';
tri=delaunay(x,y);
h = trisurf(tri, x, y, z);
l = light('Position',[-50 -15 29]);
lighting phong
shading interp
view(0, 90);
xlim([x_min x_max]);
ylim([y_min y_max]);

subplot(3,1,2)
x=[ghost.r(1,:)' ; part.r(1,:)'];
y=[ghost.r(2,:)' ; part.r(2,:)'];
z=[ghost.p' ; nan*ones(pb.N,1)];
tri=delaunay(x,y);
h = trisurf(tri, x, y, z);
%axis vis3d
l = light('Position',[-50 -15 29]);
lighting phong
shading interp
view(0, 90)
xlim([x_min x_max]);
ylim([y_min y_max]);

subplot(3,1,3)
x=[ghost.r(1,:)' ; part.r(1,:)'];
y=[ghost.r(2,:)' ; part.r(2,:)'];
z=[ghost.p' ; part.p'];
tri=delaunay(x,y);
h = trisurf(tri, x, y, z);
%axis vis3d
l = light('Position',[-50 -15 29]);
lighting phong
shading interp
view(0, 90)
xlim([x_min x_max]);
ylim([y_min y_max]);

%% Plot velocity

hf = figure;
set(hf, 'position', [220, 100, 800, 600]);

subplot(2,1,1)
x=[ghost.r(1,:)' ; part.r(1,:)'];
y=[ghost.r(2,:)' ; part.r(2,:)'];
z=[ghost.v(1,:)' ; part.v(1,:)'];
tri=delaunay(x,y);
h = trisurf(tri, x, y, z);
%axis vis3d
l = light('Position',[-50 -15 29]);
lighting phong
shading interp
view(0, 90)
xlim([x_min x_max]);
ylim([y_min y_max]);

subplot(2,1,2)
x=[ghost.r(1,:)' ; part.r(1,:)'];
y=[ghost.r(2,:)' ; part.r(2,:)'];
z=[ghost.v(2,:)' ; part.v(2,:)'];
tri=delaunay(x,y);
h = trisurf(tri, x, y, z);
%axis vis3d
l = light('Position',[-50 -15 29]);
lighting phong
shading interp
view(0, 90)
xlim([x_min x_max]);
ylim([y_min y_max]);
