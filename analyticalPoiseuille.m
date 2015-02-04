clear all
close all           
F = 5e-3;
nu = 1e-6;
L = 2e-3;
divisions = 100;
z = [0:1:divisions] / divisions;
z = z * L;

figure
for t=0:.1:2
    vx = -F / (2*nu)  * z .* (z - L);
    for n=0:20
        vx = vx - (4 * F * L ^ 2 / (nu * pi^3 * (2 * n + 1)^3)) * ...
            sin(pi * z * (2 * n + 1) / L) * ...
            exp(-(2 * n + 1) ^ 2 * pi ^ 2 * nu / L ^ 2 * t);
    end
    plot (z,vx);
    hold on;
end
hold off