function checkJacobian(obj, fd, part, edgeP, dummyP, Jf, Jg)

mCols_r = obj.cols_r;
mCols_v = obj.cols_v;
mCols_p = obj.cols_p;
pb = obj.pb;

% Finite difference Jacobians
if fd
    % perturbations
    dr = 1e-8;
    dv = 1e-5;
    dp = 1e-4;
    tic;
    [Jf_FD, Jg_FD] = jac_FD(obj, pb, part, edgeP, dummyP, [dr dv dp]);
    fprintf('Time to approximate Jacobian (FD): %f s\n', toc);
end

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
    fprintf('   2_norm entire Jacobian = %g\n', norm(Jf - Jf_FD));
    fprintf('   inf-norm entire Jacobian = %g, maxDiff = %g\n', norm(Jf - Jf_FD, inf), max(max(abs(Jf - Jf_FD))));    
    fprintf('   pos. columns    = %g\n', norm(Jf(:,mCols_r) - Jf_FD(:,mCols_r),inf));
    fprintf('   vel. columns    = %g\n', norm(Jf(:,mCols_v) - Jf_FD(:,mCols_v),inf));
    fprintf('   pres. columns   = %g\n', norm(Jf(:,mCols_p) - Jf_FD(:,mCols_p),inf));
    
    fprintf('Differences ||Jg_an - Jg_fd||\n');
    maxDiff2 = abs(Jg - Jg_FD) > .0001 * abs(Jg);
    spy(maxDiff2, 'o');
    fprintf('   entire Jacobian = %g\n', norm(Jg - Jg_FD));
    fprintf('   inf-norm entire Jacobian = %g, maxDiff = %g\n', norm(Jg - Jg_FD, inf), max(max(abs(Jg - Jg_FD))));    
    fprintf('   pos. columns    = %g\n', norm(Jg(:,mCols_r) - Jg_FD(:,mCols_r),inf));
    fprintf('   vel. columns    = %g\n', norm(Jg(:,mCols_v) - Jg_FD(:,mCols_v),inf));
    fprintf('   pres. columns   = %g\n', norm(Jg(:,mCols_p) - Jg_FD(:,mCols_p),inf));
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
g_v = Jg(:,mCols_v);
f_p = Jf(:,mCols_p);
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


