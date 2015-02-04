function [Jf, Jg] = jac(obj, pb, part, edgeP, dummyP, varargin)
%JAC  Calculates the analytical Jacobians of the index-2 Hessenberg DAE.
%     [Jf, Jg] = jac(pb, part)
%       forces a evaluation of the ghosts and of all particle neighbours.
%     [Jf, Jg] = jac(pb, part, ghosts)
%       uses the specified ghosts and assumes that the neighbours have been
%       already set (and are available in the structure 'part')

% %         if nargin == 2
% %             % Set the edgeP points.
% %             edgeP = set_ghosts(pb, part);
% % 
% %             % For each particle, find its neighbours (particles and ghosts)
% %             for i = 1 : pb.N
% %                 [nb_p, nb_d] = find_neighbours(part.r(:,i), pb, part, edgeP);
% %                 part.nb_p{i} = nb_p;
% %                 part.nb_d{i} = nb_d;
% %             end
% %         else
% %             edgeP = varargin{1};
% %         end

% Calculate the Jacobian of the RHS of the momentunm equations and the
% Jacobian of the algebraic constraints (the divergence-free conditions).
numPos = 2 * pb.N;
numVel = 2 * pb.N;
numPres = pb.N;

numVelE = 2 * pb.Ne;
numPresE = pb.Ne;

% Jacobian for particles only: rows are associated to particles only (not
% edgeP). We do not write momentum for edgeP, we write BC.
Jf = zeros(numVel,  numPos + numVel + numPres + numVelE + numPresE);
Jg = zeros(numPres, numPos + numVel + numPres + numVelE + numPresE);

for ia = 1 : pb.N
    
    % Neighbours of particle 'a'.
    nb_p = part.nb_p{a};
    nb_d = part.nb_d{a};
       
    % Rows in Jacobians corresponding to particle 'a', knowing that
    % particle 'a' is not edge
    row_Jf = ((a-1)*2 + 1 : (a-1)*2 + 2);
    row_Jg = a;

    % Columns in Jacobians corresponding to states of particle 'a'.
    col_ra = obj.Cols_r(a);
    col_va = obj.Cols_v(a);
    col_pa = obj.Cols_p(a);
    
    % Walk all particle neighbours.
    for ib = 1 : length(nb_p) + length(nb_d)
        if (ib <= length(nb_p)) % part or edgeP
            b = nb_p(ib);
            r = obj.grabR(a, part, edgeP) - obj.grabR(b, part, edgeP);
        else                    % dummyP
            dumPIdx = nb_d(ib - length(nb_p));
            b = dummyP.idx(dumPIdx); % associated edgeP particle        
        	r = obj.grabR(a, part, edgeP) - dummyP.r(:,dumPIdx);
        end

        v = obj.grabV(a, part, edgeP) - obj.grabV(b, part, edgeP);
        gradW = kernel(r, pb.h, 1);
        hessW = kernel(r, pb.h, 2);
        
        % Columns in Jacobians corresponding to particle 'b'.
        col_rb = obj.Cols_r(b);
        col_vb = obj.Cols_v(b);
        col_pb = obj.Cols_p(b);
        
        % Temporary variables.
        p_bar = obj.grabP(a, part, edgeP)/pb.rho^2 + obj.grabP(b, part, edgeP)/pb.rho^2;
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
end

return


% 
%                 
%                 
% 
% % =========================================================================
% 
% 
% % =========================================================================
% 
% function [Gr, Gv, Gp] = ghost_influence(bc)
% % The 2 dimensional array 'bc' encodes the BC for some edgeP point. We can
% % have the following cases:
% %   bc(1) = 0  or  1
% %   bc(2) = 0  or  2
% 
% if bc(2) == 0
%    Gr = [1 0; 0 1];
%    Gv = [1 0; 0 1];
% else
%    Gr = [1 0; 0 -1]; 
%    Gv = [-1 0; 0 -1];
% end
% 
% Gp = 1;
