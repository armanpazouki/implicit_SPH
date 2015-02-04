function [] = checkKernel(pb)
%CHECKKERNEL  Utility function to plot and check the 2D kernel function.

%h = 0.002;
h = pb.h;

fprintf('Using h = %10.6f\n', h);

%% Plot the kernel function
figure;
x = linspace(-2*h, 2*h, 40);
y = linspace(-2*h, 2*h, 40);
[X,Y] = meshgrid(x,y);

W = zeros(size(X));
Wx = zeros(size(X));
Wy = zeros(size(X));

for i = 1:numel(X)
    r = [X(i); Y(i)];
    W(i) = kernel(r, h, 0);
    Wr = kernel(r, h, 1);
    Wx(i) = Wr(1);
    Wy(i) = Wr(2);
end

contour(X,Y,W, 20);
hold on
quiver(X,Y,Wx,Wy);
colormap hsv
axis equal


%% Check PU property
del = h/2;
%del = pb.del;

x = (-2*h:del:2*h);
y = (-2*h:del:2*h);

% Check 2-D kernel
dA = del^2;
sum2 = 0;
count = 0;
for ix = 1:length(x)
    for iy = 1:length(y)
        r = [x(ix);y(iy)];
        w = kernel(r, h, 0);
        if w ~=0
            count = count + 1;
        end
        sum2 = sum2 + w * dA;
    end
end

fprintf('Partition of Unity check:  %f\n', sum2);
fprintf('Evaluated at %i  out of %i\n', count, length(x)*length(y));



%% Check derivatives at origin and away from it (2D).
check_derivs([1.4*h; -0.3*h], h);
check_derivs([0;0], h);

return


% =========================================================================

function [] = check_derivs(r, h)

% Analytical derivatives
g_an = kernel(r, h, 1);
H_an = kernel(r, h, 2);

% FD approximations
dx = h/100;
dy = h/100;

W = kernel(r, h, 0);

W_p1_0 = kernel(r + [dx;0], h, 0);
W_m1_0 = kernel(r + [-dx;0], h, 0);
g_forward(1) = (W_p1_0 - W) / dx;
g_central(1) = (W_p1_0 - W_m1_0) / (2 * dx);

W_0_p1 = kernel(r + [0;dy], h, 0);
W_0_m1 = kernel(r + [0;-dy], h, 0);
g_forward(2) = (W_0_p1 - W) / dy;
g_central(2) = (W_0_p1 - W_0_m1) / (2 * dy);

W_p2_0 = kernel(r + [2*dx;0], h, 0);
W_m2_0 = kernel(r + [-2*dx;0], h, 0);
H_forward(1,1) = (W_p2_0 - 2*W_p1_0 + W) / (dx*dx);
H_central(1,1) = (-W_p2_0 + 16*W_p1_0 - 30*W + 16*W_m1_0 - W_m2_0) / (12*dx^2);

W_0_p2 = kernel(r + [0;2*dy], h, 0);
W_0_m2 = kernel(r + [0;-2*dy], h, 0);
H_forward(2,2) = (W_0_p2 - 2*W_0_p1 + W) / (dy*dy);
H_central(2,2) = (-W_0_p2 + 16*W_0_p1 - 30*W + 16*W_0_m1 - W_0_m2) / (12*dy^2);

W_p1_p1 = kernel(r + [dx;dy], h, 0);
W_p1_m1 = kernel(r + [dx;-dy], h, 0);
W_m1_p1 = kernel(r + [-dx;dy], h, 0);
W_m1_m1 = kernel(r + [-dx;-dy], h, 0);
H_forward(1,2) = (W_p1_p1 - W_p1_0 - W_0_p1 + W) / (dx*dy);
H_forward(2,1) = H_forward(1,2);
H_central(1,2) = (W_p1_p1 - W_p1_m1 - W_m1_p1 + W_m1_m1) / (4*dx*dy);
H_central(2,1) = H_central(1,2);


fprintf('Derivatives for the 2-D kernel at (%f, %f)\n', r);
fprintf('   gradient\n');
fprintf('     analytical:   %10.6f  %10.6f\n', g_an(1), g_an(2));
fprintf('     forward fd:   %10.6f  %10.6f\n', g_forward(1), g_forward(2));
fprintf('     central fd:   %10.6f  %10.6f\n', g_central(1), g_central(2));
fprintf('   Hessian (xx, xy, yy)\n');
fprintf('     analytical:   %10.6f  %10.6f  %10.6f\n', H_an(1,1), H_an(1,2), H_an(2,2));
fprintf('     forward fd:   %10.6f  %10.6f  %10.6f\n', H_forward(1,1), H_forward(1,2), H_forward(2,2));
fprintf('     central fd:   %10.6f  %10.6f  %10.6f\n', H_central(1,1), H_central(1,2), H_central(2,2));



