%CHECKGHOSTS  Check and plot the particles and ghosts for a 2-D problem.
%   The check consists of calculating the SPH density estimate at various
%   points and comparing it against the base value.

% Get the problem dimension from the user:
while true
    dim_str = input('Dimension (1, 2, or 3): ', 's');
    dim = str2double(dim_str);
    if (~isempty(dim) && isnumeric(dim) && (dim == 1 || dim == 2 || dim == 3))
        break;
    end
end

pb = init_problem(dim);
part = init_particles(pb);

% Set the ghost points.
ghost = set_ghosts(pb, part);

% Plot particles and ghosts
hf = figure;
set(hf, 'position', [50, 50, 640, 640]);
hold on

plot_particles(hf, part, ghost, pb);

% Find and plot neighbors for a few locations.
% Estimate density at those locations.
fprintf('Specified density:  %f\n', pb.rho);

switch pb.dim
    case 1
        % center location.
        loc = 0;
        rho = calc_density(loc, pb, part, ghost);
        fprintf('  at %f  rho = %f\n', loc, rho);
        
        % location close to the right boundary.
        loc = 2.5*pb.h;
        rho = calc_density(loc, pb, part, ghost);
        fprintf('  at %f  rho = %f\n', loc, rho);
        
        % particle close to left boundary.
        loc = part.r(2);
        plot_neighbours(hf, loc, pb, part, ghost);
        plot(loc, 0, 'r*', 'markerfacecolor', 'r');
        rho = calc_density(loc, pb, part, ghost);
        fprintf('  at %f  rho = %f\n', loc, rho);
        
    case 2
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
        ip = round(pb.n(1)/2);
        jp = pb.n(2)-1;
        loc = part.r(:,(ip-1)*pb.n(2) + jp);
        plot_neighbours(hf, loc, pb, part, ghost);
        plot(loc(1), loc(2), 'r*', 'markerfacecolor', 'r');
        rho = calc_density(loc, pb, part, ghost);
        fprintf('  at [%f %f]  rho = %f\n', loc(1), loc(2), rho);
        
    case 3
        % center location.
        loc = [0;0;0];
        rho = calc_density(loc, pb, part, ghost);
        fprintf('  at [%f %f %f]  rho = %f\n', loc(1), loc(2), loc(3), rho);

        % particle close to left boundary.
        loc = part.r(:,5);
        plot_neighbours(hf, loc, pb, part, ghost);
        plot3(loc(1), loc(2), loc(3), 'r*', 'markerfacecolor', 'r');
        rho = calc_density(loc, pb, part, ghost);
        fprintf('  at [%f %f %f]  rho = %f\n', loc(1), loc(2), loc(3), rho);
end
