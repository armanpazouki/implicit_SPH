function [Jf, Jg] = jac(pb, part, varargin)
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
    for i = 1 : part.num
        [nb_p, nb_g] = find_neighbours(part.r(:,i), pb, part, ghost);
        part.nb_p{i} = nb_p;
        part.nb_g{i} = nb_g;
    end
else
    ghost = varargin{1};
end

% Calculate the Jacobian of the RHS of the momentunm equations and the
% Jacobian of the algebraic constraints (the divergence-free conditions).
numPos = pb.dim * part.num;
numVel = pb.dim * part.num;
numPres = part.num;

Jf = zeros(numVel,  numPos + numVel + numPres);
Jg = zeros(numPres, numPos + numVel + numPres);

for a = 1 : part.num
    
    % Neighbours of particle 'a'.
    nb_p = part.nb_p{a};
    nb_g = part.nb_g{a};
       
    % Rows in Jacobians corresponding to particle 'a'.
    row_Jf = ((a-1)*pb.dim + 1 : (a-1)*pb.dim + pb.dim);
    row_Jg = a;

    % Columns in Jacobians corresponding to states of particle 'a'.
    col_ra = cols_r(a, pb.dim, part.num);
    col_va = cols_v(a, pb.dim, part.num);
    col_pa = cols_p(a, pb.dim, part.num);
    
    % Walk all particle neighbours.
    for ib = 1 : length(nb_p)
        b = nb_p(ib);
        r = part.r(:,a) - part.r(:,b);
        v = part.v(:,a) - part.v(:,b);
        gradW = kernel(r, pb.h, 1);
        hessW = kernel(r, pb.h, 2);
        
        % Columns in Jacobians corresponding to particle 'b'.
        col_rb = cols_r(b, pb.dim, part.num);
        col_vb = cols_v(b, pb.dim, part.num);
        col_pb = cols_p(b, pb.dim, part.num);
        
        % Temporary variables.
        p_bar = part.p(a)/part.rho(a)^2 + part.p(b)/part.rho(b)^2;
        rho_bar = (part.rho(a) + part.rho(b)) / 2;
        mu_bar = part.mu(a) + part.mu(b);
        rr_bar = r' * r + pb.eta2;
       
        % Jacobian blocks corresponding to the term A.
        Aa_ra = pb.m * p_bar * hessW;
        Jf(row_Jf, col_ra) = Jf(row_Jf, col_ra) - Aa_ra;
        Jf(row_Jf, col_rb) = Jf(row_Jf, col_rb) + Aa_ra;
        
        Aa_pa = pb.m / part.rho(a)^2 * gradW';
        Aa_pb = pb.m / part.rho(b)^2 * gradW';
        Jf(row_Jf, col_pa) = Jf(row_Jf, col_pa) - Aa_pa;
        Jf(row_Jf, col_pb) = Jf(row_Jf, col_pb) - Aa_pb;
        
        % Jacobian blocks corresponding to the term B.
        Ba_ra = pb.m * (mu_bar / rho_bar^2 / rr_bar) * ...
            (v * r' * (hessW - 2 * (r' * gradW') / rr_bar * eye(pb.dim)) + v * gradW);
        Jf(row_Jf, col_ra) = Jf(row_Jf, col_ra) + Ba_ra;
        Jf(row_Jf, col_rb) = Jf(row_Jf, col_rb) - Ba_ra;
        
        Ba_va = pb.m * (mu_bar / rho_bar^2 / rr_bar) * (r' * gradW') * eye(pb.dim);
        Jf(row_Jf, col_va) = Jf(row_Jf, col_va) + Ba_va;
        Jf(row_Jf, col_vb) = Jf(row_Jf, col_vb) - Ba_va;
        
        % Jacobian blocks corresponding to the term C.
        Ca_ra = pb.m / part.rho(b) * v' * hessW;
        Jg(row_Jg, col_ra) = Jg(row_Jg, col_ra) + Ca_ra;
        Jg(row_Jg, col_rb) = Jg(row_Jg, col_rb) - Ca_ra;
        
        Ca_va = pb.m / part.rho(b) * gradW;
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
        col_rb = cols_r(ghost.idx(b), pb.dim, part.num);
        col_vb = cols_v(ghost.idx(b), pb.dim, part.num);
        col_pb = cols_p(ghost.idx(b), pb.dim, part.num);
        
        % Temporary variables.
        p_bar = part.p(a)/part.rho(a)^2 + ghost.p(b)/ghost.rho(b)^2;
        rho_bar = (part.rho(a) + ghost.rho(b)) / 2;
        mu_bar = part.mu(a) + ghost.mu(b);
        rr_bar = r' * r + pb.eta2;
  
        % Get relationship ghost - associated particle
        [Gr, Gv] = ghost_influence(pb.dim, ghost.bc(:,b));
        
        % Jacobian blocks corresponding to the term A.
        Aa_ra = pb.m * p_bar * hessW;
        Jf(row_Jf, col_ra) = Jf(row_Jf, col_ra) - Aa_ra;
        Jf(row_Jf, col_rb) = Jf(row_Jf, col_rb) + Aa_ra * Gr;
        
        Aa_pa = pb.m / part.rho(a)^2 * gradW';
        Aa_pb = pb.m / ghost.rho(b)^2 * gradW';
        Jf(row_Jf, col_pa) = Jf(row_Jf, col_pa) - Aa_pa;
        Jf(row_Jf, col_pb) = Jf(row_Jf, col_pb) - Aa_pb;
        
        % Jacobian blocks corresponding to the term B.
        Ba_ra = pb.m * (mu_bar / rho_bar^2 / rr_bar) * ...
            (v * r' * (hessW - 2 * (r' * gradW') / rr_bar * eye(pb.dim)) + v * gradW);
        Jf(row_Jf, col_ra) = Jf(row_Jf, col_ra) + Ba_ra;
        Jf(row_Jf, col_rb) = Jf(row_Jf, col_rb) - Ba_ra * Gr;
        
        Ba_va = pb.m * (mu_bar / rho_bar^2 / rr_bar) * (r' * gradW') * eye(pb.dim);
        Jf(row_Jf, col_va) = Jf(row_Jf, col_va) + Ba_va;
        Jf(row_Jf, col_vb) = Jf(row_Jf, col_vb) - Ba_va * Gv;

        % Jacobian blocks corresponding to the term C.
        Ca_ra = pb.m / ghost.rho(b) * v' * hessW;
        Jg(row_Jg, col_ra) = Jg(row_Jg, col_ra) + Ca_ra;
        Jg(row_Jg, col_rb) = Jg(row_Jg, col_rb) - Ca_ra * Gr;
        
        Ca_va = pb.m / ghost.rho(b) * gradW;
        Jg(row_Jg, col_va) = Jg(row_Jg, col_va) + Ca_va;
        Jg(row_Jg, col_vb) = Jg(row_Jg, col_vb) - Ca_va * Gv;
    end
    
    
end

return

% =========================================================================

function c = cols_r(i, dim, n)
c = ((i-1)*dim + 1 : i*dim);
return

function c = cols_v(i, dim, n)
numPos = dim * n;
c = (numPos + (i-1)*dim + 1 : numPos + i*dim);
return

function c = cols_p(i, dim, n)
c = 2 * dim * n + i;
return

% =========================================================================

function [Gr, Gv] = ghost_influence(dim, bc)

r = ones(dim, 1);
ii = find(bc == 2);
r(ii) = -1;
Gr = diag(r);

if isempty(ii)
    Gv = eye(dim);
else
    Gv = zeros(dim);
end

