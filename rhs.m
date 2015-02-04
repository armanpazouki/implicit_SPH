function [f, g] = rhs(pb, part, ghost, varargin)
%RHS  Evaluate the right hand side terms in the index-2 Hessenberg DAE.
%     [f, g] = rhs(pb, part)
%        forces a evaluaion of the ghosts and of all particle neighbours.
%     [f, g] = rhs(pb, part, ghosts)
%       uses the specified ghosts and assumes that the neighbours have been
%       already set (and are available in the structure 'part').

    %     if nargin == 2
    %         % Set the ghost points.
    %         ghost = set_ghosts(pb, part);
    % 
    %         % For each particle, find its neighbours (particles and ghosts)
    %         for i = 1 : pb.N
    %             [nb_p, nb_g] = find_neighbours(part.r(:,i), pb, part, ghost);
    %             part.nb_p{i} = nb_p;
    %             part.nb_g{i} = nb_g;
    %         end
    %     else
    %         ghost = varargin{1};
    %     end

% Evaluate the RHS of the momentum equations and the divergence-free
% algebraic constraints.
f = zeros(2, pb.N);
g = zeros(1, pb.N);

for a = 1 : pb.N
    A = zeros(2, 1);
    B = zeros(2, 1);
    C = 0;
   
%     fprintf('\n-----------\na = %i\n', a);
    
    nb_p = part.nb_p{a};
    nb_g = part.nb_g{a};
    
    for ib = 1 : length(nb_p)
        b = nb_p(ib);
        r = part.r(:,a) - part.r(:,b);
        v = part.v(:,a) - part.v(:,b);
        gradW = kernel(r, pb.h, 1);
        
        p_bar = part.p(a)/pb.rho^2 + part.p(b)/pb.rho^2;
        rho_bar = (pb.rho + pb.rho) / 2;
        mu_bar = pb.mu + pb.mu;
        rr_bar = r' * r + pb.eta2;
      
        Aa = pb.m * p_bar * gradW';
        Ba = pb.m * (mu_bar / rho_bar^2 / rr_bar) * (r' * gradW') * v;
        Ca = pb.m / pb.rho * gradW * v;
        
        A = A + Aa;
        B = B + Ba;
        C = C + Ca;
        
%         fprintf('   np: %3i   q=%10.6f  gradW=[%10.6f  %10.6f]  r=[%10.6f  %10.6f]  v=[%10.6f  %10.6f]', b, norm(r)/pb.h, gradW, r, v);
%         fprintf('   (%10.6f %10.6f %10.6f %10.6f )', p_bar, rho_bar, mu_bar, rr_bar);
%         fprintf('   [%10.6f  %10.6f]  [%10.6f  %10.6f]\n', Aa, Ba);
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
        
        Aa = pb.m * p_bar * gradW';
        Ba = pb.m * (mu_bar / rho_bar^2 / rr_bar) * (r' * gradW') * v;
        Ca = pb.m / pb.rho * gradW * v;
        
        A = A + Aa;
        B = B + Ba;
        C = C + Ca;
        
%         fprintf('   ng: %3i   q=%10.6f  gradW=[%10.6f  %10.6f]  r=[%10.6f  %10.6f]  v=[%10.6f  %10.6f]', ghost.idx(b), norm(r)/pb.h, gradW, r, v);
%         fprintf('   (%10.6f %10.6f %10.6f %10.6f )', p_bar, rho_bar, mu_bar, rr_bar);
%         fprintf('   [%10.6f  %10.6f]  [%10.6f  %10.6f]\n', Aa, Ba);
    end
    
%     fprintf('A = %g  %g\n', A);
%     fprintf('B = %g  %g\n', B);
    
    f(:, a) = [pb.F;0] - A + B;

    g(a) = C;
end


% Reshape 'f' and 'g' as column vectors.
f = reshape(f, 2 * pb.N, 1);
g = g';

