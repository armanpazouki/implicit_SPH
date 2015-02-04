%CHECKJACOBIAN  Check analytical Jacobian against finite difference approx.
%   checkJacobian calculates the analytical Jacobian for a 2-D problem and
%   estimates it using forward finite differences.

%% Evaluate FD Jacobian?
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

%% Initialize problem and particles
pb = init_problem();
part = init_particles(pb);

%% Set the ghost points and the neighbours for all particles
% (this is not really needed here, but we want to have them for various
% other checks).
ghost = set_ghosts(pb, part);

for i = 1 : pb.N
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

numPos = 2 * pb.N;
numVel = 2 * pb.N;
numPres = pb.N;

cols_r = (1:numPos);
cols_v = (numPos+1:numPos+numVel);
cols_p = (numPos+numVel+1:numPos+numVel+numPres);

% Difference in Jacobian of momentum RHS
fprintf('\n');
fprintf('Dimensions\n');
fprintf('   num particles     = %i\n', pb.N);
fprintf('   part. in each dim = %i %i\n', pb.nx, pb.ny);
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

    maxDiff1 = abs(Jf - Jf_FD) > (2000000 * abs(Jf));
    spy(abs(Jf_FD) > .000001, 'o');
    max(max(Jf))
    fprintf('   entire Jacobian = %g\n', norm(Jf - Jf_FD));
%     fprintf('   entire Jacobian = %g, maxDiff = %g\n', norm(Jf - Jf_FD), norm(Jf - Jf_FD, inf));
    fprintf('   pos. columns    = %g\n', norm(Jf(:,cols_r) - Jf_FD(:,cols_r)));
    fprintf('   vel. columns    = %g\n', norm(Jf(:,cols_v) - Jf_FD(:,cols_v)));
    fprintf('   pres. columns   = %g\n', norm(Jf(:,cols_p) - Jf_FD(:,cols_p)));
    
    fprintf('Differences ||Jg_an - Jg_fd||\n');
    maxDiff2 = abs(Jg - Jg_FD) > .0001 * abs(Jg);
    spy(maxDiff2, 'o');
    fprintf('   entire Jacobian = %g\n', norm(Jg - Jg_FD));
    fprintf('   pos. columns    = %g\n', norm(Jg(:,cols_r) - Jg_FD(:,cols_r)));
    fprintf('   vel. columns    = %g\n', norm(Jg(:,cols_v) - Jg_FD(:,cols_v)));
    fprintf('   pres. columns   = %g\n', norm(Jg(:,cols_p) - Jg_FD(:,cols_p)));
    fprintf('\n');
    
    %% max error
    relativeErrorJf = zeros(size(Jf));
    nonZeroIndNJf = find(Jf > 1e-10);
    relative_error_non_zoro_componentsJf = (Jf(nonZeroIndNJf) - Jf_FD(nonZeroIndNJf)) ./ Jf(nonZeroIndNJf);
    zeroIndNJf = find(Jf==0);
    supposedlyZeroComponentsJf = Jf(zeroIndNJf) - Jf_FD(zeroIndNJf);
    fprintf('   max Jf = %g\n', max(max(Jf)));
    fprintf('   max relative errorJf 1 = %g\n', max(max(relative_error_non_zoro_componentsJf)));
    fprintf('   max errorJf 1 = %g\n', max(max(supposedlyZeroComponentsJf)));
    
    relativeErrorJg = zeros(size(Jg));
    nonZeroIndNJg = find(Jg > 1e-10);
    relative_error_non_zoro_componentsJg = (Jg(nonZeroIndNJg) - Jg_FD(nonZeroIndNJg)) ./ Jg(nonZeroIndNJg);
    zeroIndNJg = find(Jg==0);
    supposedlyZeroComponentsJg = Jg(zeroIndNJg) - Jg_FD(zeroIndNJg);
    fprintf('   max Jg = %g\n', max(max(Jg)));
    fprintf('   max relative errorJg 1 = %g\n', max(max(relative_error_non_zoro_componentsJg)));
    fprintf('   max errorJg 1 = %g\n', max(max(supposedlyZeroComponentsJg)));
    %%
end

fprintf('Rank of Jg\n');
rr = rank(Jg);
fprintf('   rank(Jg_an) = %i    deficiency = %i\n', rr, pb.N - rr);

if fd
    rr_FD = rank(Jg_FD);
    fprintf('   rank(Jg_fd) = %i    deficiency = %i\n', rr_FD, pb.N - rr_FD);
end

fprintf('\n');
fprintf('Hessenberg index-2 condition\n');
g_v = Jg(:,cols_v);
f_p = Jf(:,cols_p);
rr = rank(g_v * f_p);
passed = (rr == pb.N);
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
set(gca, 'xtick', (2.5 : 2 : numVel-2.5),'xticklabel', []);
set(gca, 'ytick', (1.5 : 1 : numPres-0.5), 'yticklabel',[]);
set(gca,'xcolor',[0.7 0.7 0.7], 'ycolor', [0.7 0.7 0.7]);
xlabel('');
grid on
title('g_v')

% Split g_v by components of the velocities.
figure

for i = 1 : 2
    subplot(2, 1, i)
    g_vi = g_v(:,i:2:numVel);
    spy(g_vi)
    set(gca, 'xticklabel', [], 'yticklabel',[]);
    xlabel('');
    tstr = sprintf('g_{v%i}', i);
    title(tstr);
end

