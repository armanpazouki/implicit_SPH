function ghost = set_ghosts(pb, part)
x_min = 0;
x_max = pb.Lx;
Dx = x_max - x_min;

y_min = -(0.5*pb.Ly);
y_max = (0.5*pb.Ly);
Dy = y_max - y_min;

%% BC in x direction

% "Minus" boundary
idx_pm = find(part.r(1,:) > x_max - 2 * pb.h);
ng = length(idx_pm);
r_pm = [2 * x_max - part.r(1, idx_pm) ; part.r(2, idx_pm)];  % wall
v_pm = repmat(2*[pb.uR;pb.vR], 1, ng) - part.v(:,idx_pm);
p_pm = part.p(idx_pm);
bc_pm = repmat([2;0], 1, ng);
% r_pm = part.r(:, idx_pm) + repmat([-Dx;0], 1, ng); % periodic
% v_pm = part.v(:,idx_pm);
% p_pm = part.p(idx_pm) + pb.K * pb.Lx;
% bc_pm = repmat([-1;0], 1, ng);

% "Plus" boundary
idx_pp = find(part.r(1,:) < x_min + 2 * pb.h);
ng = length(idx_pp);
r_pp = [2 * x_min - part.r(1, idx_pp) ; part.r(2, idx_pp)];  % wall
v_pp = repmat(2*[pb.uL;pb.vL], 1, ng) - part.v(:,idx_pp);
p_pp = part.p(idx_pp);
bc_pp = repmat([-2;0], 1, ng);
% r_pp = part.r(:, idx_pp) + repmat([Dx;0], 1, ng);
% v_pp = part.v(:,idx_pp);
% p_pp = part.p(idx_pp) - pb.K * pb.Lx;
% bc_pp = repmat([1;0], 1, ng);

% Ghosts for x direction BC
gx.idx = [idx_pm idx_pp];
gx.r = [r_pm r_pp];
gx.v = [v_pm v_pp];
gx.p = [p_pm p_pp];
gx.bc = [bc_pm bc_pp];

%% BC in y direction

% "Minus" boundary
idx_pm = find(part.r(2,:) < y_min + 2 * pb.h);
ng = length(idx_pm);
r_pm = [part.r(1, idx_pm) ; 2 * y_min - part.r(2, idx_pm)];
% v_pm = zeros(2, ng);
v_pm = repmat(2*[pb.uB;pb.vB], 1, ng) - part.v(:, idx_pm);
p_pm = part.p(idx_pm);
bc_pm = repmat([0;-2], 1, ng);

idx_gm = find(gx.r(2,:) < y_min + 2 * pb.h);
ng = length(idx_gm);
r_gm = [gx.r(1, idx_gm) ; 2 * y_min - gx.r(2, idx_gm)];
v_gm = zeros(2, ng);
% v_gm = repmat(2*[pb.uB;pb.vB], 1, ng) - gx.v(:, idx_gm);
p_gm = gx.p(idx_gm);
% bc_gm = repmat([1;2], 1, ng);
bc_gm = repmat([2;2], 1, ng);

% "Plus" boundary
idx_pp = find(part.r(2,:) > y_max - 2 * pb.h);
ng = length(idx_pp);
r_pp = [part.r(1, idx_pp) ; 2 * y_max - part.r(2, idx_pp)];
% v_pp = zeros(2, ng);
v_pp = repmat(2*[pb.uT;pb.vT], 1, ng) - part.v(:, idx_pp);
p_pp = part.p(idx_pp);
bc_pp = repmat([0;2], 1, ng);

idx_gp = find(gx.r(2,:) > y_max - 2 * pb.h);
ng = length(idx_gp);
r_gp = [gx.r(1, idx_gp) ; 2 * y_max - gx.r(2, idx_gp)];
v_gp = zeros(2, ng);
% v_gp = repmat(2*[pb.uT;pb.vT], 1, ng) - gx.v(:, idx_gp);
p_gp = gx.p(:, idx_gp);
% bc_gp = repmat([1;2], 1, ng);
bc_gp = repmat([2;2], 1, ng);

% Ghosts for y direction BC
gy.idx = [idx_pm gx.idx(idx_gm) idx_pp gx.idx(idx_gp)];
gy.r = [r_pm r_gm r_pp r_gp];
gy.v = [v_pm v_gm v_pp v_gp];
gy.p = [p_pm p_gm p_pp p_gp];
gy.bc = [bc_pm bc_gm bc_pp bc_gp];


%% Collect ghost data

ghost.idx = [gx.idx gy.idx];
ghost.r = [gx.r gy.r];
ghost.v = [gx.v gy.v];
ghost.p = [gx.p gy.p];
ghost.bc = [gx.bc gy.bc];

ghost.num = length(ghost.idx);


% ghost.vv = zeros(2,ghost.num);
% K = pb.K + pb.rho * pb.F;
% tmp = K * (0.5*pb.Ly)^2 / (2 * pb.mu);
% 
% for i = 1:ghost.num
%     ghost.vv(1,i) = tmp * (1-(ghost.r(2,i)/(0.5*pb.Ly))^2);
% end



return
