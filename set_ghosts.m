function [edgeP, dummyP] = set_ghosts(pb, part)
x_min = 0;
x_max = pb.Lx;
Dx = x_max - x_min;

y_min = 0;%-(0.5*pb.Ly);
y_max = pb.Ly;
Dy = y_max - y_min;

%% bottom edge and dummy
%edge
x = pb.del * (0:pb.nx+1);
sizeEdge1 = size(x,2);
y = 0 * ones(1, sizeEdge1);

idx_ed1 = [1:sizeEdge1];
ed1_r = [x;y];
ed1_v = [pb.uB * ones(1, sizeEdge1); pb.vB * ones(1, sizeEdge1)];
ed1_p = 0 * ones(1, sizeEdge1);
bc_e1 = repmat([0;-2], 1, sizeEdge1);

%dummy
yDum = [-3:-1] * pb.del;
[X,Y] = meshgrid(x,yDum);
numDum1 = sizeEdge1 * size(yDum, 2);
X = reshape(X, 1, numDum1);
Y = reshape(Y, 1, numDum1);

dum1_r = [X;Y];
dum1_v = zeros(2, size(X,2));
dum1_p = zeros(1, size(X,2));
[idx_dum1,Y] = meshgrid(idx_ed1,yDum);
idx_dum1 = reshape(idx_dum1, 1, numDum1);
bc_d1 = repmat([0;-3], 1, numDum1);

%% top edge
%edge
x = pb.del * (0:pb.nx+1);
sizeEdge2 = size(x,2);
y = (pb.ny + 1) * pb.del * ones(1, sizeEdge2);

idx_ed2 = [1:sizeEdge2] + sizeEdge1;
ed2_r = [x;y];
ed2_v = [pb.uT * ones(1, sizeEdge2); pb.vT * ones(1, sizeEdge2)];
ed2_p = 0 * ones(1, sizeEdge2);
bc_e2 = repmat([0;2], 1, sizeEdge2);

%dummy
yDum = [pb.ny + 2:pb.ny + 4] * pb.del;
[X,Y] = meshgrid(x,yDum);
numDum2 = sizeEdge2 * size(yDum, 2);
X = reshape(X, 1, numDum2);
Y = reshape(Y, 1, numDum2);

dum2_r = [X;Y];
dum2_v = zeros(2, size(X,2));
dum2_p = zeros(1, size(X,2));
[idx_dum2,Y] = meshgrid(idx_ed2,yDum);
idx_dum2 = reshape(idx_dum2, 1, numDum2);
bc_d2 = repmat([0;3], 1, numDum2);

%% left edge
%edge
y = pb.del * (1:pb.ny);
sizeEdge3 = size(y,2);
x = 0 * ones(1, sizeEdge3);

idx_ed3 = [1:sizeEdge3] + sizeEdge1 + sizeEdge2;
ed3_r = [x;y];
ed3_v = [pb.uL * ones(1, sizeEdge3); pb.vL * ones(1, sizeEdge3)];
ed3_p = 0 * ones(1, sizeEdge3);
bc_e3 = repmat([-2;0], 1, sizeEdge3);

%dummy
xDum = [-3:-1] * pb.del;
[X,Y] = meshgrid(xDum,y);
numDum3 = sizeEdge3 * size(xDum, 2);
X = reshape(X, 1, numDum3);
Y = reshape(Y, 1, numDum3);

dum3_r = [X;Y];
dum3_v = zeros(2, size(X,2));
dum3_p = zeros(1, size(X,2));
[X,idx_dum3] = meshgrid(xDum,idx_ed3);
idx_dum3 = reshape(idx_dum3, 1, numDum3);
bc_d3 = repmat([-3;0], 1, numDum3);

%% right edge
%edge
y = pb.del * (1:pb.ny);
sizeEdge4 = size(y,2);
x = (pb.nx + 1) * pb.del * ones(1, sizeEdge4);

idx_ed4 = [1:sizeEdge4] + sizeEdge1 + sizeEdge2 + sizeEdge3;
ed4_r = [x;y];
ed4_v = [pb.uR * ones(1, sizeEdge4); pb.vR * ones(1, sizeEdge4)];
ed4_p = 0 * ones(1, sizeEdge4);
bc_e4 = repmat([2;0], 1, sizeEdge4);

%dummy
xDum = [pb.nx+2:pb.nx+4] * pb.del;
[X,Y] = meshgrid(xDum,y);
numDum4 = sizeEdge4 * size(xDum, 2);
X = reshape(X, 1, numDum4);
Y = reshape(Y, 1, numDum4);

dum4_r = [X;Y];
dum4_v = zeros(2, size(X,2));
dum4_p = zeros(1, size(X,2));
[X,idx_dum4] = meshgrid(xDum,idx_ed4);
idx_dum4 = reshape(idx_dum4, 1, numDum4);
bc_d4 = repmat([3;0], 1, numDum4);

% % bc type
% bc_pm = repmat([2;0], 1, ng);
% gx.bc = [bc_pm bc_pp];

%% corner dummies
% TR-->5
xDum = [pb.nx+2:pb.nx+4] * pb.del;
yDum = [pb.ny+1:pb.ny+4] * pb.del;
[X,Y] = meshgrid(xDum,yDum);
numDum5 = size(X, 1) * size(X, 2);
X = reshape(X, 1, numDum5);
Y = reshape(Y, 1, numDum5);

dum5_r = [X;Y];
dum5_v = zeros(2, size(X,2));
dum5_p = zeros(1, size(X,2));
idx_dum5 = idx_ed2(sizeEdge2) * ones(size(X)); % same index for all, equal to TR edge particle
bc_d5 = repmat([3;0], 1, numDum5);

% TL-->6
xDum = [-3:-1] * pb.del;
yDum = [pb.ny+1:pb.ny+4] * pb.del;
[X,Y] = meshgrid(xDum,yDum);
numDum6 = size(X, 1) * size(X, 2);
X = reshape(X, 1, numDum6);
Y = reshape(Y, 1, numDum6);

dum6_r = [X;Y];
dum6_v = zeros(2, size(X,2));
dum6_p = zeros(1, size(X,2));
idx_dum6 = idx_ed2(1) * ones(size(X)); % same index for all, equal to TR edge particle
bc_d6 = repmat([3;0], 1, numDum6);

% BL-->7
xDum = [-3:-1] * pb.del;
yDum = [-3:0] * pb.del;
[X,Y] = meshgrid(xDum,yDum);
numDum7 = size(X, 1) * size(X, 2);
X = reshape(X, 1, numDum7);
Y = reshape(Y, 1, numDum7);

dum7_r = [X;Y];
dum7_v = zeros(2, size(X,2));
dum7_p = zeros(1, size(X,2));
idx_dum7 = idx_ed1(1) * ones(size(X)); % same index for all, equal to TR edge particle
bc_d7 = repmat([3;0], 1, numDum7);

% BR-->8
xDum = [pb.nx+2:pb.nx+4] * pb.del;
yDum = [-3:0] * pb.del;
[X,Y] = meshgrid(xDum,yDum);
numDum8 = size(X, 1) * size(X, 2);
X = reshape(X, 1, numDum8);
Y = reshape(Y, 1, numDum8);

dum8_r = [X;Y];
dum8_v = zeros(2, size(X,2));
dum8_p = zeros(1, size(X,2));
idx_dum8 = idx_ed1(sizeEdge1) * ones(size(X)); % same index for all, equal to TR edge particle
bc_d8 = repmat([3;0], 1, numDum8);

%% Collect edge data

edgeP.idx = [idx_ed1 idx_ed2 idx_ed3 idx_ed4];
edgeP.r = [ed1_r ed2_r ed3_r ed4_r];
edgeP.v = [ed1_v ed2_v ed3_v ed4_v];
edgeP.p = [ed1_p ed2_p ed3_p ed4_p];
edgeP.bc = [bc_e1 bc_e2 bc_e3 bc_e4];

edgeP.num = length(edgeP.idx);

% edgeP.nb_p = cell(edgeP.num, 1);
% edgeP.nb_e = cell(edgeP.num, 1);
% edgeP.nb_d = cell(edgeP.num, 1);

%% Collect dummy data

dummyP.idx = [idx_dum1 idx_dum2 idx_dum3 idx_dum4 idx_dum5 idx_dum6 idx_dum7 idx_dum8];
dummyP.r = [dum1_r dum2_r dum3_r dum4_r dum5_r dum6_r dum7_r dum8_r];
dummyP.v = [dum1_v dum2_v dum3_v dum4_v dum5_v dum6_v dum7_v dum8_v];
dummyP.p = [dum1_p dum2_p dum3_p dum4_p dum5_p dum6_p dum7_p dum8_p];
dummyP.bc = [bc_d1 bc_d2 bc_d3 bc_d4 bc_d5 bc_d6 bc_d7 bc_d8];

dummyP.num = length(dummyP.idx);





return
