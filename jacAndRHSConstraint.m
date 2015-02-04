function [Jc, c] = jacAndRHSConstraint(obj, pb, part, ghost, varargin)

% Calculate the Jacobian of the constraints as well as the right hand side of the momentunm equations and the
% Jacobian of the algebraic constraints (the divergence-free conditions).

    numPos = obj.numPos;
    numVel = obj.numVel;
    numPres = obj.numPres;

    numGhost = ghost.num;
    numVelG = obj.numVel_ghost;
    numPresG = obj.numPres_ghost;

    Jc = zeros(numVelG + numPresG + 1,  numPos + numVel + numPres + numVelG + numPresG);
    c = zeros(numVelG + numPresG + 1, 1);
    
    numSavedVelConstraint = 2 * numGhost;
    bcType = [0;0];
    % Walk all ghost neighbours.
    for i = 1 : numGhost
        bcType = ghost.bc(:,i);
        origPartInd = ghost.idx(i);
        
        v_a = ghost.v(:,i); 
        v_b = part.v(:,origPartInd);
        p_a = ghost.p(i); 
        p_b = part.p(origPartInd);
        
        col_va = obj.Cols_v_ghost(i, pb.N, numGhost);
        col_pa = obj.Cols_p_ghost(i, pb.N, numGhost);
        col_vb = obj.Cols_v(origPartInd, pb.N);
        col_pb = obj.Cols_p(origPartInd, pb.N);
        
        row_Jc_v = (2*(i-1)+1:2*(i-1)+2);
        row_Jc_p = numSavedVelConstraint + i;
        
        if isequal(bcType,[1;0]) % periodic, ghosts at the inlet
            Jc(row_Jc_v, col_va) = eye(2);
            Jc(row_Jc_v, col_vb) = -eye(2);
            Jc(row_Jc_p, col_pa) = 1;
            Jc(row_Jc_p, col_pb) = -1;

            c(row_Jc_v) = v_a - v_b;
            c(row_Jc_p) = p_a - p_b - pb.K * pb.Lx;
            
        elseif isequal(bcType,[-1;0]) % periodic, ghosts at the outlet
            Jc(row_Jc_v, col_va) = eye(2);
            Jc(row_Jc_v, col_vb) = -eye(2);
            Jc(row_Jc_p, col_pa) = 1;
            Jc(row_Jc_p, col_pb) = -1;

            c(row_Jc_v) = v_a - v_b;
            c(row_Jc_p) = p_a - p_b + pb.K * pb.Lx;
            
        elseif isequal(bcType,[0;2]) % top is wall
            Jc(row_Jc_v, col_va) = eye(2);
            Jc(row_Jc_v, col_vb) = eye(2);
            Jc(row_Jc_p, col_pa) = 1;
            Jc(row_Jc_p, col_pb) = -1;

            c(row_Jc_v) = v_a + v_b - 2 * [pb.uT; pb.vT];
            c(row_Jc_p) = p_a - p_b;
            
        elseif isequal(bcType,[0;-2]) % bottom is wall
            Jc(row_Jc_v, col_va) = eye(2);
            Jc(row_Jc_v, col_vb) = eye(2);
            Jc(row_Jc_p, col_pa) = 1;
            Jc(row_Jc_p, col_pb) = -1;

            c(row_Jc_v) = v_a + v_b - 2 * [pb.uB; pb.vB];
            c(row_Jc_p) = p_a - p_b;
            
        elseif isequal(bcType,[2;0]) % Right is wall
            Jc(row_Jc_v, col_va) = eye(2);
            Jc(row_Jc_v, col_vb) = eye(2);
            Jc(row_Jc_p, col_pa) = 1;
            Jc(row_Jc_p, col_pb) = -1;

            c(row_Jc_v) = v_a + v_b - 2 * [pb.uR; pb.vR];
            c(row_Jc_p) = p_a - p_b;
            
        elseif isequal(bcType,[-2;0])
            Jc(row_Jc_v, col_va) = eye(2);
            Jc(row_Jc_v, col_vb) = eye(2);
            Jc(row_Jc_p, col_pa) = 1;
            Jc(row_Jc_p, col_pb) = -1;

            c(row_Jc_v) = v_a + v_b - 2 * [pb.uL; pb.vL];
            c(row_Jc_p) = p_a - p_b;
            
        elseif isequal(bcType,[1;2])
            Jc(row_Jc_v, col_va) = eye(2);
            Jc(row_Jc_v, col_vb) = eye(2);
            Jc(row_Jc_p, col_pa) = 1;
            Jc(row_Jc_p, col_pb) = -1;

            c(row_Jc_v) = v_a + v_b;
            c(row_Jc_p) = p_a - p_b;
            
        elseif isequal(bcType,[2;2])
            Jc(row_Jc_v, col_va) = eye(2);
            Jc(row_Jc_p, col_pa) = 1;
            Jc(row_Jc_p, col_pb) = -1;

            c(row_Jc_v) = v_a;
            c(row_Jc_p) = p_a - p_b;
        
        else
            'not defined boundary condition'
        end
    end
    row_Jc_p = numSavedVelConstraint + numGhost + 1;
    col_p = obj.Cols_p(1, pb.N);
    Jc(row_Jc_p, col_p) = 1;
    c(row_Jc_p) = 0;
return