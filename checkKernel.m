function [] = checkKernel()
%CHECKKERNEL  Utility function to plot anf check the kernel function.
%    It plots the kernel and its first derivative for 1D, 2D, and 3D, and
%    performs some sanity checks on the kernel function.

h = 1;

fprintf('Using h = %10.6f\n', h);

% Plot the kernel functions for 1D, 2D, 3D
hf = figure;
set(hf, 'position', [50, 50, 1000, 300]);

subplot(1,3,1)
plotKernel_1D(h);
subplot(1,3,2)
plotKernel_2D(h);
subplot(1,3,3)
plotKernel_3D(h);

% Check PU for 1D, 2D, 3D
check_PU(h);

% Check derivatives at origin and away from it (2D).
check_derivs([1.4*h; -0.3*h], h);
check_derivs([0;0], h);

return

% =========================================================================

function [] = plotKernel_1D(h)

x = linspace(-2*h, 2*h, 40);
W = zeros(size(x));
Wx = zeros(size(x));
for i = 1:numel(x)
    r = x(i);
    W(i) = kernel(r, h, 0);
    Wx(i) = kernel(r, h, 1);
end
hold on
plot(x, W, 'r')
plot(x, Wx, 'b');
grid on
box on

return

% =========================================================================

function [] = plotKernel_2D(h)

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

return

% =========================================================================

function [] = plotKernel_3D(h)

x = linspace(-2*h, 2*h, 20);
y = linspace(-2*h, 2*h, 20);
z = linspace(-2*h, 2*h, 20);
[X,Y,Z] = meshgrid(x,y,z);

Wx = zeros(size(X));
Wy = zeros(size(X));
Wz = zeros(size(X));

for i = 1:numel(X)
    r = [X(i); Y(i); Z(i)];
    Wr = kernel(r, h, 1);
    Wx(i) = Wr(1);
    Wy(i) = Wr(2);
    Wz(i) = Wr(3);
end

quiver3(X,Y,Z,Wx,Wy,Wz);

% =========================================================================

function [] = check_PU(h)

del = h/10;
x = (-2*h:del:2*h);
y = (-2*h:del:2*h);
z = (-2*h:del:2*h);

% Check 1-D kernel
dx = del;
sum1 = 0;
for ix = 1:length(x)
   r = x(ix);
   w = kernel(r, h, 0);
   sum1 = sum1 + w * dx;
end

% Check 2-D kernel
dA = del^2;
sum2 = 0;
for ix = 1:length(x)
    for iy = 1:length(y)
        r = [x(ix);y(iy)];
        w = kernel(r, h, 0);
        sum2 = sum2 + w * dA;
    end
end

% Check 3-D kernel
dV = del^3;
sum3 = 0;
for ix = 1:length(x)
    for iy = 1:length(y)
        for iz = 1:length(z)
            r = [x(ix);y(iy);z(iz)];
            w = kernel(r, h, 0);
            sum3 = sum3 + w * dV;
        end
    end
end

fprintf('Partition of Unity check\n');
fprintf('   1-D:  %f\n', sum1);
fprintf('   2-D:  %f\n', sum2);
fprintf('   3-D:  %f\n', sum3);

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



