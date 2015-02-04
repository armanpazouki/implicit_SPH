function [Jf, Jg] = jac_FD(pb, part, deltas)
%JAC_FD  Finite-diff approximation to the Jacobians of the index-2 Hessenberg DAE.

% perturbations
dr = deltas(1);
dv = deltas(2);
dp = deltas(3);

numPos = 2 * pb.N;
numVel = 2 * pb.N;
numPres = pb.N;

Jf = zeros(numVel,  numPos + numVel + numPres);
Jg = zeros(numPres, numPos + numVel + numPres);

% Evaluate the FD Jacobian approximations, one column at a time.
% Note that we force a ghost re-evaluation at every RHS evaluation (to
% ensure accuracy of the FD approximation).
[f,g] = rhs(pb, part);

for a = 1 : pb.N
    % Derivatives w.r.t. r(a)
    for i = 1 : 2
        col = (a-1) * 2 + i;
        part.r(i,a) = part.r(i,a) + dr;
        [ff, gg] = rhs(pb, part);
        Jf(:,col) = (ff - f) / dr;
        Jg(:,col) = (gg - g) / dr;
        part.r(i,a) = part.r(i,a) - dr;
    end
    
    % Derivatives w.r.t. v(a)
    for i = 1 : 2
        col = numPos + (a-1) * 2 + i;
        part.v(i,a) = part.v(i,a) + dv;
        [ff, gg] = rhs(pb, part);
        Jf(:,col) = (ff - f) / dv;
        Jg(:,col) = (gg - g) / dv;
        part.v(i,a) = part.v(i,a) - dv;
    end
    
    % Derivatives w.r.t p(a)
    col = numPos + numVel + a;
    part.p(a) = part.p(a) + dp;
    [ff, gg] = rhs(pb, part);
    Jf(:,col) = (ff - f) / dp;
    Jg(:,col) = (gg - g) / dp;
    part.p(a) = part.p(a) - dp;
end

