function [Jf, Jg] = jacB(pb, part, varargin)
%JAC  Calculates the analytical Jacobians of the index-2 Hessenberg DAE.
%     [Jf, Jg] = jac(pb, part)
%       forces a evaluation of the ghosts and of all particle neighbours.
%     [Jf, Jg] = jac(pb, part, ghosts)
%       uses the specified ghosts and assumes that the neighbours have been
%       already set (and are available in the structure 'part')

if nargin == 2
    % Set the ghost points.
    ghost = set_ghosts(pb, part);
    
    % For each particle, find its neighbours (particles and ghosts)
    for i = 1 : pb.N
        [nb_p, nb_g] = find_neighbours(part.r(:,i), pb, part, ghost);
        part.nb_p{i} = nb_p;
        part.nb_g{i} = nb_g;
    end
else
    ghost = varargin{1};
end

% Calculate the Jacobian of the RHS of the momentunm equations and the
% Jacobian of the algebraic constraints (the divergence-free conditions).
numPos = 2 * pb.N;
numVel = 2 * pb.N;
numPres = pb.N;

Jf = zeros(numVel,  numPos + numVel + numPres);
Jg = zeros(numPres, numPos + numVel + numPres);

for a = 1 : pb.N
    
    % Neighbours of particle 'a'.
    nb_p = part.nb_p{a};
    nb_g = part.nb_g{a};
       
    % Rows in Jacobians corresponding to particle 'a'.
    row_Jf = ((a-1)*2 + 1 : (a-1)*2 + 2);
    row_Jg = a;

    % Columns in Jacobians corresponding to states of particle 'a'.
    col_ra = cols_r(a, pb.N);
    col_va = cols_v(a, pb.N);
    col_pa = cols_p(a, pb.N);
    
    % Walk all particle neighbours.
    for ib = 1 : length(nb_p)
        b = nb_p(ib);
        r = part.r(:,a) - part.r(:,b);
        v = part.v(:,a) - part.v(:,b);
        gradW = kernel(r, pb.h, 1);
        hessW = kernel(r, pb.h, 2);
        
        % Columns in Jacobians corresponding to particle 'b'.
        col_rb = cols_r(b, pb.N);
        col_vb = cols_v(b, pb.N);
        col_pb = cols_p(b, pb.N);
        
        % Temporary variables.
        p_bar = part.p(a)/pb.rho^2 + part.p(b)/pb.rho^2;
        rho_bar = (pb.rho + pb.rho) / 2;
        mu_bar = pb.mu + pb.mu;
        rr_bar = r' * r + pb.eta2;
       
        % Jacobian blocks corresponding to the term A.
        Aa_ra = pb.m * p_bar * hessW;
        Jf(row_Jf, col_ra) = Jf(row_Jf, col_ra) - Aa_ra;
        Jf(row_Jf, col_rb) = Jf(row_Jf, col_rb) + Aa_ra;
        
        Aa_pa = pb.m / pb.rho^2 * gradW';
        Aa_pb = pb.m / pb.rho^2 * gradW';
        Jf(row_Jf, col_pa) = Jf(row_Jf, col_pa) - Aa_pa;
        Jf(row_Jf, col_pb) = Jf(row_Jf, col_pb) - Aa_pb;
        
        % Jacobian blocks corresponding to the term B.
        Ba_ra = pb.m * (mu_bar / rho_bar^2 / rr_bar) * ...
            (v * r' * (hessW - 2 * (r' * gradW') / rr_bar * eye(2)) + v * gradW);
        Jf(row_Jf, col_ra) = Jf(row_Jf, col_ra) + Ba_ra;
        Jf(row_Jf, col_rb) = Jf(row_Jf, col_rb) - Ba_ra;
        
        Ba_va = pb.m * (mu_bar / rho_bar^2 / rr_bar) * (r' * gradW') * eye(2);
        Jf(row_Jf, col_va) = Jf(row_Jf, col_va) + Ba_va;
        Jf(row_Jf, col_vb) = Jf(row_Jf, col_vb) - Ba_va;
        
        % Jacobian blocks corresponding to the term C.
        Ca_ra = pb.m / pb.rho * v' * hessW;
        Jg(row_Jg, col_ra) = Jg(row_Jg, col_ra) + Ca_ra;
        Jg(row_Jg, col_rb) = Jg(row_Jg, col_rb) - Ca_ra;
        
        Ca_va = pb.m / pb.rho * gradW;
        Jg(row_Jg, col_va) = Jg(row_Jg, col_va) + Ca_va;
        Jg(row_Jg, col_vb) = Jg(row_Jg, col_vb) - Ca_va;
    end
    
    % Walk all ghost neighbours.
    for ib = 1 : length(nb_g)
        b = nb_g(ib);
        r = part.r(:,a) - ghost.r(:,b);
        v = part.v(:,a) - ghost.v(:,b);
        gradW = kernel(r, pb.h, 1);
        hessW = kernel(r, pb.h, 2);
        
        % Columns in Jacobians corresponding to particle ghost.idx(b).
        col_rb = cols_r(ghost.idx(b), pb.N);
        col_vb = cols_v(ghost.idx(b), pb.N);
        col_pb = cols_p(ghost.idx(b), pb.N);
        
        % Temporary variables.
        p_bar = part.p(a)/pb.rho^2 + ghost.p(b)/pb.rho^2;
        rho_bar = (pb.rho + pb.rho) / 2;
        mu_bar = pb.mu + pb.mu;
        rr_bar = r' * r + pb.eta2;
  
        % Get relationship ghost - associated particle
        [Gr, Gv, Gp] = ghost_influence(ghost.bc(:,b));
        
        % Jacobian blocks corresponding to the term A.
        Aa_ra = pb.m * p_bar * hessW;
        Jf(row_Jf, col_ra) = Jf(row_Jf, col_ra) - Aa_ra;
        Jf(row_Jf, col_rb) = Jf(row_Jf, col_rb) + Aa_ra * Gr;
        
        Aa_pa = pb.m / pb.rho^2 * gradW';
        Aa_pb = pb.m / pb.rho^2 * gradW';
        Jf(row_Jf, col_pa) = Jf(row_Jf, col_pa) - Aa_pa;
        Jf(row_Jf, col_pb) = Jf(row_Jf, col_pb) - Aa_pb * Gp;
        
        % Jacobian blocks corresponding to the term B.
        Ba_ra = pb.m * (mu_bar / rho_bar^2 / rr_bar) * ...
            (v * r' * (hessW - 2 * (r' * gradW') / rr_bar * eye(2)) + v * gradW);
        Jf(row_Jf, col_ra) = Jf(row_Jf, col_ra) + Ba_ra;
        Jf(row_Jf, col_rb) = Jf(row_Jf, col_rb) - Ba_ra * Gr;
        
        Ba_va = pb.m * (mu_bar / rho_bar^2 / rr_bar) * (r' * gradW') * eye(2);
        Jf(row_Jf, col_va) = Jf(row_Jf, col_va) + Ba_va;
        Jf(row_Jf, col_vb) = Jf(row_Jf, col_vb) - Ba_va * Gv;

        % Jacobian blocks corresponding to the term C.
        Ca_ra = pb.m / pb.rho * v' * hessW;
        Jg(row_Jg, col_ra) = Jg(row_Jg, col_ra) + Ca_ra;
        Jg(row_Jg, col_rb) = Jg(row_Jg, col_rb) - Ca_ra * Gr;
        
        Ca_va = pb.m / pb.rho * gradW;
        Jg(row_Jg, col_va) = Jg(row_Jg, col_va) + Ca_va;
        Jg(row_Jg, col_vb) = Jg(row_Jg, col_vb) - Ca_va * Gv;
    end
    
    
end

return

% =========================================================================

function c = cols_r(i, n)
c = ((i-1)*2 + 1 : i*2);
return

function c = cols_v(i, n)
numPos = 2 * n;
c = (numPos + (i-1)*2 + 1 : numPos + i*2);
return

function c = cols_p(i, n)
c = 2 * 2 * n + i;
return

% =========================================================================

function [Gr, Gv, Gp] = ghost_influence(bc)
% The 2 dimensional array 'bc' encodes the BC for some ghost point. We can
% have the following cases:
%   bc(1) = 0  or  1
%   bc(2) = 0  or  2

if bc(2) == 0
   Gr = [1 0; 0 1];
   Gv = [1 0; 0 1];
else
   Gr = [1 0; 0 -1]; 
   Gv = [-1 0; 0 -1];
end

Gp = 1;
