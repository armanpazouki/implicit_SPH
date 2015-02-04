function ghost = set_ghosts(pb, part)
%
% NOTE:  Currently, we assume periodic BC all around!!!
%

switch pb.dim
    case 1
        ghost = set_ghosts_1D(pb, part);
    case 2
        ghost = set_ghosts_2D(pb, part);
    case 3
        ghost = set_ghosts_3D(pb, part);
end


return 

% =========================================================================

function ghost = set_ghosts_1D(pb, part)
%
% Always using periodic BC in x.
%

x_min = pb.domain(1,1);
x_max = pb.domain(1,2);
Dx = x_max - x_min;

% Region -x
idx_m = find(part.r(1,:) > x_max - 2 * pb.h);
ng = length(idx_m);
r_m = part.r(:, idx_m) + repmat(-Dx, 1, ng);
v_m = part.v(:, idx_m);
p_m = part.p(idx_m);
bc_m = ones(1, ng);

% Region +x
idx_p = find(part.r(1,:) < x_min + 2 * pb.h);
ng = length(idx_p);
r_p = part.r(:, idx_p) + repmat(Dx, 1, ng);
v_p = part.v(:, idx_p);
p_p = part.p(idx_p);
bc_p = ones(1, ng);

% Collect information in the structure 'g'.
ghost.idx = [idx_m  idx_p];
ghost.r = [r_m  r_p];
ghost.v = [v_m  v_p];
ghost.p = [p_m  p_p];
ghost.bc = [bc_m bc_p];

ghost.num = length(ghost.idx);

ghost.rho = part.rho(ghost.idx);
ghost.mu = part.mu(ghost.idx);

return

% =========================================================================

function ghost = set_ghosts_2D(pb, part)

x_min = pb.domain(1,1);
x_max = pb.domain(1,2);
Dx = x_max - x_min;

y_min = pb.domain(2,1);
y_max = pb.domain(2,2);
Dy = y_max - y_min;

% BC in x direction.
switch pb.BC(1)

    case 1  % Periodic
        
        % "Minus" boundary
        idx_pm = find(part.r(1,:) > x_max - 2 * pb.h);
        ng = length(idx_pm);
        r_pm = part.r(:, idx_pm) + repmat([-Dx;0], 1, ng);
        v_pm = part.v(:,idx_pm);
        p_pm = part.p(idx_pm);
        bc_pm = repmat([1;0], 1, ng);
        
        % "Plus" boundary
        idx_pp = find(part.r(1,:) < x_min + 2 * pb.h);
        ng = length(idx_pp);
        r_pp = part.r(:, idx_pp) + repmat([Dx;0], 1, ng);
        v_pp = part.v(:,idx_pp);
        p_pp = part.p(idx_pp);
        bc_pp = repmat([1;0], 1, ng);
        
        % Ghosts for x direction BC
        gx.idx = [idx_pm idx_pp];
        gx.r = [r_pm r_pp];
        gx.v = [v_pm v_pp];
        gx.p = [p_pm p_pp];
        gx.bc = [bc_pm bc_pp];
        
    case 2  % Wall
        
        % "Minus" boundary
        idx_pm = find(part.r(1,:) < x_min + 2 * pb.h);
        ng = length(idx_pm);
        r_pm = [2 * x_min - part.r(1, idx_pm) ; part.r(2, idx_pm)];
        v_pm = zeros(2, ng);
        p_pm = part.p(idx_pm);
        bc_pm = repmat([2;0], 1, ng);
               
        % "Plus" boundary
        idx_pp = find(part.r(1,:) > x_max - 2 * pb.h);
        ng = length(idx_pp);
        r_pp = [2 * x_max - part.r(1, idx_pp) ; part.r(2, idx_pp)];
        v_pp = zeros(2, ng);
        p_pp = part.p(idx_pp);
        bc_pp = repmat([2;0], 1, ng);

        % Ghosts for x direction BC
        gx.idx = [idx_pm idx_pp];
        gx.r = [r_pm r_pp];
        gx.v = [v_pm v_pp];
        gx.p = [p_pm p_pp];
        gx.bc = [bc_pm bc_pp];        
        
end

% BC in y direction
switch pb.BC(2)
    
    case 1  % Periodic
        
        % "Minus" boundary
        idx_pm = find(part.r(2,:) > y_max - 2 * pb.h);
        ng = length(idx_pm);
        r_pm = part.r(:, idx_pm) + repmat([0;-Dy], 1, ng);
        v_pm = part.v(:,idx_pm);
        p_pm = part.p(idx_pm);
        bc_pm = repmat([0;1], 1, ng);
        
        idx_gm = find(gx.r(2,:) > y_max - 2 * pb.h);
        ng = length(idx_gm);
        r_gm = gx.r(:, idx_gm) + repmat([0;-Dy], 1, ng);
        v_gm = gx.v(:, idx_gm);
        p_gm = gx.p(idx_gm);
        bc_gm = repmat([pb.BC(1);1], 1, ng);
        
        % "Plus" boundary
        idx_pp = find(part.r(2,:) < y_min + 2 * pb.h);
        ng = length(idx_pp);
        r_pp = part.r(:, idx_pp) + repmat([0;Dy], 1, ng);
        v_pp = part.v(:,idx_pp);
        p_pp = part.p(idx_pp);
        bc_pp = repmat([0;1], 1, ng);
        
        idx_gp = find(gx.r(2,:) < y_min + 2 * pb.h);
        ng = length(idx_gp);
        r_gp = gx.r(:, idx_gp) + repmat([0;Dy], 1, ng);
        v_gp = gx.v(:, idx_gp);
        p_gp = gx.p(:, idx_gp);
        bc_gp = repmat([pb.BC(1);1], 1, ng);
        
        % Ghosts for y direction BC
        gy.idx = [idx_pm gx.idx(idx_gm) idx_pp gx.idx(idx_gp)];
        gy.r = [r_pm r_gm r_pp r_gp];
        gy.v = [v_pm v_gm v_pp v_gp];
        gy.p = [p_pm p_gm p_pp p_gp];
        gy.bc = [bc_pm bc_gm bc_pp bc_gp];
        
    case 2  % Wall
        
        % "Minus" boundary
        idx_pm = find(part.r(2,:) < y_min + 2 * pb.h);
        ng = length(idx_pm);
        r_pm = [part.r(1, idx_pm) ; 2 * y_min - part.r(2, idx_pm)];
        v_pm = zeros(2, ng);
        p_pm = part.p(idx_pm);
        bc_pm = repmat([0;2], 1, ng);
       
        idx_gm = find(gx.r(2,:) < y_min + 2 * pb.h);
        ng = length(idx_gm);
        r_gm = [gx.r(1, idx_gm) ; 2 * y_min - gx.r(2, idx_gm)];
        v_gm = zeros(2, ng);
        p_gm = gx.p(idx_gm);
        bc_gm = repmat([pb.BC(1);2], 1, ng);
        
        % "Plus" boundary
        idx_pp = find(part.r(2,:) > y_max - 2 * pb.h);
        ng = length(idx_pp);
        r_pp = [part.r(1, idx_pp) ; 2 * y_max - part.r(2, idx_pp)];
        v_pp = zeros(2, ng);
        p_pp = part.p(idx_pp);
        bc_pp = repmat([0;2], 1, ng);
        
        idx_gp = find(gx.r(2,:) > y_max - 2 * pb.h);
        ng = length(idx_gp);
        r_gp = [gx.r(1, idx_gp) ; 2 * y_max - gx.r(2, idx_gp)];
        v_gp = zeros(2, ng);
        p_gp = gx.p(:, idx_gp);
        bc_gp = repmat([pb.BC(1);2], 1, ng);

        % Ghosts for y direction BC
        gy.idx = [idx_pm gx.idx(idx_gm) idx_pp gx.idx(idx_gp)];
        gy.r = [r_pm r_gm r_pp r_gp];
        gy.v = [v_pm v_gm v_pp v_gp];
        gy.p = [p_pm p_gm p_pp p_gp];
        gy.bc = [bc_pm bc_gm bc_pp bc_gp];
      
end

ghost.idx = [gx.idx gy.idx];
ghost.r = [gx.r gy.r];
ghost.v = [gx.v gy.v];
ghost.p = [gx.p gy.p];
ghost.bc = [gx.bc gy.bc];

ghost.num = length(ghost.idx);

ghost.rho = part.rho(ghost.idx);
ghost.mu = part.mu(ghost.idx);


return

% =========================================================================

function ghost = set_ghosts_3D(pb, part)

% TODO
% Currently only periodic BC in all dircetions!!!
%

x_min = pb.domain(1,1);
x_max = pb.domain(1,2);
Dx = x_max - x_min;

y_min = pb.domain(2,1);
y_max = pb.domain(2,2);
Dy = y_max - y_min;

z_min = pb.domain(3,1);
z_max = pb.domain(3,2);
Dz = z_max - z_min;

ng = 0;

for i = 1 : part.num
    r = part.r(:,i);
    g = zeros(3,1);
    dupX = false;
    dupY = false;
    dupZ = false;
    
    if r(1) > x_max - 2*pb.h
        dupX = true;
        g(1) = r(1) - Dx;
    end
    if r(1) < x_min + 2*pb.h
        dupX = true;
        g(1) = r(1) + Dx;
    end
    
    if r(2) > y_max - 2*pb.h
        dupY = true;
        g(2) = r(2) - Dy;
    end
    if r(2) < y_min + 2*pb.h
        dupY = true;
        g(2) = r(2) + Dy;
    end
    
    if r(3) > z_max - 2*pb.h
        dupZ = true;
        g(3) = r(3) - Dz;
    end
    if r(3) < z_min + 2*pb.h
        dupZ = true;
        g(3) = r(3) + Dz;
    end
    
    if dupX
        % ghost in x,0,0
        ng = ng + 1;
        ghost.r(:,ng) = [g(1);r(2);r(3)];
        ghost.idx(ng) = i;
        
        if dupY
            % ghost in x,y,0
            ng = ng + 1;
            ghost.r(:,ng) = [g(1);g(2);r(3)];
            ghost.idx(ng) = i;
            
            if dupZ
                % ghost in x,y,z
                ng = ng + 1;
                ghost.r(:,ng) = [g(1);g(2);g(3)];
                ghost.idx(ng) = i;
                
            end
        end
        
        if dupZ
            % ghost in x,0,z
            ng = ng + 1;
            ghost.r(:,ng) = [g(1);r(2);g(3)];
            ghost.idx(ng) = i;
            
        end
    end
    
    if dupY
        % ghost in 0,y,0
        ng = ng + 1;
        ghost.r(:,ng) = [r(1);g(2);r(3)];
        ghost.idx(ng) = i;
        
        if dupZ
            % ghost in 0,y,z
            ng = ng + 1;
            ghost.r(:,ng) = [r(1);g(2);g(3)];
            ghost.idx(ng) = i;
            
        end
    end
    
    if dupZ
        % ghost in 0,0,z
        ng = ng + 1;
        ghost.r(:,ng) = [r(1);r(2);g(3)];
        ghost.idx(ng) = i;
        
    end
    

end % end for

ghost.num = length(ghost.idx);

ghost.rho = part.rho(ghost.idx);
ghost.mu = part.mu(ghost.idx);
ghost.p = part.p(ghost.idx);
ghost.v = part.v(:,ghost.idx);

% TODO
ghost.bc = repmat([1;1;1], 1, ghost.num);

return


