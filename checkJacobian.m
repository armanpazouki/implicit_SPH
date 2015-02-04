%CHECKJACOBIAN  Check analytical Jacobian against finite difference approx.
%   checkJacobian calculates the analytical Jacobian for a 2-D problem and
%   estimates it using forward finite differences.

% Get the problem dimension from the user:
while true
    dim_str = input('Dimension (1, 2, or 3): ', 's');
    dim = str2double(dim_str);
    if (~isempty(dim) && isnumeric(dim) && (dim == 1 || dim == 2 || dim == 3))
        break;
    end
end

% Evaluate FD Jacobian?
fd = true;
while true
   fd_str = input('Evaluate FD Jacobian? (Y/N): ', 's');
   fd_str = upper(fd_str);
   if (~isempty(fd_str) && fd_str(1) == 'Y')
       fd = true;
       break;
   end
   if (~isempty(fd_str) && fd_str(1) == 'N')
       fd = false;
       break;
   end
end


pb = init_problem(dim);
part = init_particles(pb);

% Set the ghost points and the neighbours for all particles
% (this is not really needed here, but we want to have them for various
% other checks).
ghost = set_ghosts(pb, part);

for i = 1 : part.num
    [nb_p, nb_g] = find_neighbours(part.r(:,i), pb, part, ghost);
    part.nb_p{i} = nb_p;
    part.nb_g{i} = nb_g;
end

% Base value of the RHS (force ghosts re-evaluation)
tic;
[f,g] = rhs(pb, part);
fprintf('Time to evaluate function:         %f s\n', toc);

% Analytical Jacobian (force ghosts re-evaluation)
tic;
[Jf, Jg] = jac(pb, part);
fprintf('Time to evaluate Jacobian:         %f s\n', toc);

% Finite difference Jacobians
if fd
    % perturbations
    dr = 1e-8;
    dv = 1e-5;
    dp = 1e-4;
    tic;
    [Jf_FD, Jg_FD] = jac_FD(pb, part, [dr dv dp]);
    fprintf('Time to approximate Jacobian (FD): %f s\n', toc);
end

numPos = pb.dim * part.num;
numVel = pb.dim * part.num;
numPres = part.num;

cols_r = (1:numPos);
cols_v = (numPos+1:numPos+numVel);
cols_p = (numPos+numVel+1:numPos+numVel+numPres);

% Difference in Jacobian of momentum RHS
fprintf('\n');
fprintf('Dimensions\n');
fprintf('   dim               = %i-D\n', pb.dim);
fprintf('   num particles     = %i\n', part.num);
fprintf('   part. in each dim =');
fprintf(' %i', pb.n);
fprintf('\n');
BCtypes = {'periodic', 'wall'};
fprintf('   BC in each dim    =');
fprintf(' %s', BCtypes{pb.BC});
fprintf('\n');
fprintf('   size of Jf        = %ix%i\n', size(Jf));
fprintf('   size of Jg        = %ix%i\n', size(Jg));
fprintf('\n');

if fd
    fprintf('Perturbation parameters\n');
    fprintf('   position = %e\n', dr);
    fprintf('   velocity = %e\n', dv);
    fprintf('   pressure = %e\n', dp);
    fprintf('\n');
    fprintf('Differences ||Jf_an - Jf_fd||\n');
    fprintf('   entire Jacobian = %g\n', norm(Jf - Jf_FD));
    fprintf('   pos. columns    = %g\n', norm(Jf(:,cols_r) - Jf_FD(:,cols_r)));
    fprintf('   vel. columns    = %g\n', norm(Jf(:,cols_v) - Jf_FD(:,cols_v)));
    fprintf('   pres. columns   = %g\n', norm(Jf(:,cols_p) - Jf_FD(:,cols_p)));
    fprintf('Differences ||Jg_an - Jg_fd||\n');
    fprintf('   entire Jacobian = %g\n', norm(Jg - Jg_FD));
    fprintf('   pos. columns    = %g\n', norm(Jg(:,cols_r) - Jg_FD(:,cols_r)));
    fprintf('   vel. columns    = %g\n', norm(Jg(:,cols_v) - Jg_FD(:,cols_v)));
    fprintf('   pres. columns   = %g\n', norm(Jg(:,cols_p) - Jg_FD(:,cols_p)));
    fprintf('\n');
end

fprintf('Rank of Jg\n');
rr = rank(Jg);
fprintf('   rank(Jg_an) = %i    deficiency = %i\n', rr, part.num - rr);

if fd
    rr_FD = rank(Jg_FD);
    fprintf('   rank(Jg_fd) = %i    deficiency = %i\n', rr_FD, part.num - rr_FD);
end

fprintf('\n');
fprintf('Hessenberg index-2 condition\n');
g_v = Jg(:,cols_v);
f_p = Jf(:,cols_p);
rr = rank(g_v * f_p);
passed = (rr == part.num);
fprintf('   rank(g_v * f_p) = %i   pass? %i\n', rr, passed);
fprintf('   rank(g_v) = %i\n', rank(g_v));
fprintf('   rank(f_p) = %i\n', rank(f_p));
fprintf('\n');

% Plot sparsity pattern of the sub-Jacobian g_v (i.e. the Jacobian of
% the algebraic constraints w.r.t. velocities).
figure
spy(g_v)
% spy(g_v, 's', 10);
% markerH = findall(gca,'color','b');
% set(markerH,'MarkerFaceColor','r');
set(gca, 'GridLineStyle', '-');
set(gca, 'xlim',[0.5 numVel+0.5], 'ylim',[0.5 numPres+0.5]);
set(gca, 'xtick', (pb.dim+0.5 : pb.dim : numVel-pb.dim+0.5),'xticklabel', []);
set(gca, 'ytick', (1.5 : 1 : numPres-0.5), 'yticklabel',[]);
set(gca,'xcolor',[0.7 0.7 0.7], 'ycolor', [0.7 0.7 0.7]);
xlabel('');
grid on
title('g_v')

% Split g_v by components of the velocities.
if pb.dim > 1
    figure
    
    for i = 1 : pb.dim
        subplot(pb.dim, 1, i)
        g_vi = g_v(:,i:2:numVel);
        spy(g_vi)
        set(gca, 'xticklabel', [], 'yticklabel',[]);
        xlabel('');
        tstr = sprintf('g_{v%i}', i);
        title(tstr);
    end
end

