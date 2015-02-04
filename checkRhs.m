%CHECKRHS  Check RHS at a point in the middle of the domain.

%% Initialize problem, particles, ghosts
pb = init_problem();
part = init_particles(pb);
ghost = set_ghosts(pb, part);

%% Set the neighbours for all particles
for i = 1 : pb.N
    [nb_p, nb_g] = find_neighbours(part.r(:,i), pb, part, ghost);
    part.nb_p{i} = nb_p;
    part.nb_g{i} = nb_g;
end

%% Select a point in the middle of the domain.
i = floor((pb.nx+1)/2);
j = floor((pb.ny+1)/2);
%i = 1;
%j = 1;
a = (i-1) * pb.ny + j;

nb_p = part.nb_p{a};
nb_g = part.nb_g{a};

fprintf('Evaluating RHS for particle (%i,%i)\n', i, j);
fprintf('Number of particle neighbors: %i\n', numel(nb_p));
fprintf('Number of ghost neighbors: %i\n', numel(nb_g));

A = zeros(2, 1);
B = zeros(2, 1);
C = 0;

%     fprintf('\n-----------\na = %i\n', a);

for ib = 1 : length(nb_p)
    b = nb_p(ib);
    r = part.r(:,a) - part.r(:,b);
    v = part.v(:,a) - part.v(:,b);
    gradW = kernel(r, pb.h, 1);
        
    Aa = (part.p(a) + part.p(b)) * gradW';
    Ba = (gradW * r) / (r' * r + pb.eta2) * v;
    Ca = gradW * v;
    
    A = A + Aa;
    B = B + Ba;
    C = C + Ca;
end

for ib = 1 : length(nb_g)
    b = nb_g(ib);
    r = part.r(:,a) - ghost.r(:,b);
    v = part.v(:,a) - ghost.v(:,b);
    gradW = kernel(r, pb.h, 1);
    
    p_bar = part.p(a)/pb.rho^2 + ghost.p(b)/pb.rho^2;
    rho_bar = (pb.rho + pb.rho) / 2;
    mu_bar = pb.mu + pb.mu;
    rr_bar = r' * r + pb.eta2;
    
    Aa = (part.p(a) + ghost.p(b)) * gradW';
    Ba = (gradW * r) / (r' * r + pb.eta2) * v;
    Ca = gradW * v;
    
    A = A + Aa;
    B = B + Ba;
    C = C + Ca;
end

A = (pb.m / pb.rho^2) * A;
B = (pb.m * 2 * pb.mu / pb.rho^2) * B;
C = (pb.m / pb.rho) * C;

f = B - A + [pb.F ; 0];
g = C;

