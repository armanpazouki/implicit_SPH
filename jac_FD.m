function [Jf, Jg] = jac_FD(obj, pb, part, edgeP, dummyP, deltas)
%JAC_FD  Finite-diff approximation to the Jacobians of the index-2 Hessenberg DAE.

% perturbations
dr = deltas(1);
dv = deltas(2);
dp = deltas(3);

numPos = 2 * pb.N;
numVel = 2 * pb.N;
numPres = pb.N;
%%
% Jacobian for particles only: rows are associated to particles only (not
% edgeP). We do not write momentum for edgeP, we write BC.
Jf = zeros(numVel,  numPos + numVel + numPres);
Jg = zeros(numPres, numPos + numVel + numPres);

[f,g] = rhs(obj, pb, part, edgeP, dummyP);

for a = 1 : pb.N
    % Columns in Jacobians corresponding to states of particle 'a'.
    col_ra = obj.Cols_r(a);
    col_va = obj.Cols_v(a);
    col_pa = obj.Cols_p(a);
    
    % Derivatives w.r.t. r(a)
    for i = 1 : 2
        col = col_ra(i);
        part.r(i,a) = part.r(i,a) + dr;
        [ff, gg] = rhs(obj, pb, part, edgeP, dummyP);
        Jf(:,col) = (ff - f) / dr;
        Jg(:,col) = (gg - g) / dr;
        part.r(i,a) = part.r(i,a) - dr;
    end
    
    % Derivatives w.r.t. v(a)
    for i = 1 : 2
        col = col_va(i);
        part.v(i,a) = part.v(i,a) + dv;
        [ff, gg] = rhs(obj, pb, part, edgeP, dummyP);
        Jf(:,col) = (ff - f) / dv;
        Jg(:,col) = (gg - g) / dv;
        part.v(i,a) = part.v(i,a) - dv;
    end
    
    % Derivatives w.r.t p(a)
    col = col_pa;
    part.p(a) = part.p(a) + dp;
    [ff, gg] = rhs(obj, pb, part, edgeP, dummyP);
    Jf(:,col) = (ff - f) / dp;
    Jg(:,col) = (gg - g) / dp;
    part.p(a) = part.p(a) - dp;
end
