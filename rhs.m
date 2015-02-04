function [f, g] = rhs(obj, pb, part, edgeP, dummyP)
%RHS  Evaluate the right hand side terms in the index-2 Hessenberg DAE.
%     [f, g] = rhs(pb, part)
%        forces a evaluaion of the ghosts and of all particle neighbours.
%     [f, g] = rhs(pb, part, ghosts)
%       uses the specified ghosts and assumes that the neighbours have been
%       already set (and are available in the structure 'part').


% Evaluate the RHS of the momentum equations and the divergence-free
% algebraic constraints.
f = zeros(2, pb.N);
g = zeros(1, pb.N);

for a = 1 : pb.N
    A = zeros(2, 1);
    B = zeros(2, 1);
    C = 0;
   
%     fprintf('\n-----------\na = %i\n', a);
    
    nb_p = obj.nb_p{a};
    nb_d = obj.nb_d{a};
    
    counterV = 0;
    for ib = 1 : length(nb_p) + length(nb_d)
        if (ib <= length(nb_p)) % part or edgeP
            dumPIdx = -1;
            b = nb_p(ib);
            if (b == a)
%                 continue;
            end
            r = obj.grabR(a, part, edgeP) - obj.grabR(b, part, edgeP);
        else                    % dummyP
            dumPIdx = nb_d(ib - length(nb_p));
            b = dummyP.idx(dumPIdx); % associated edgeP particle 
            b = b + pb.N;            % to comply with the grabV and grabP functions
        	r = obj.grabR(a, part, edgeP) - dummyP.r(:,dumPIdx);
        end
        
        v = obj.grabV(a, part, edgeP) - obj.grabV(b, part, edgeP);
        gradW = kernel(r, pb.h, 1);
        
        p_bar = obj.grabP(a, part, edgeP)/pb.rho^2 + obj.grabP(b, part, edgeP)/pb.rho^2;
        rho_bar = (pb.rho + pb.rho) / 2;
        mu_bar = pb.mu + pb.mu;
        rr_bar = r' * r + pb.eta2;
      
        Aa = pb.m * p_bar * gradW';
        Ba = pb.m * (mu_bar / rho_bar^2 / rr_bar) * (r' * gradW') * v;
        Ca = pb.m / pb.rho * gradW * v;
        
        A = A + Aa;
        B = B + Ba;
        C = C + Ca;
%         D = D + pb.m / pb.rho * kernel(r, pb.h, 0) * v;
        
%         vv = obj.grabV(b, part, edgeP);
%         if (vv(1) ~= 0 | vv(2) ~= 0)
%             counterV = counterV + 1;
%         end
        
%         fprintf('   np: %3i   q=%10.6f  gradW=[%10.6f  %10.6f]  r=[%10.6f  %10.6f]  v=[%10.6f  %10.6f]', b, norm(r)/pb.h, gradW, r, v);
%         fprintf('   (%10.6f %10.6f %10.6f %10.6f )', p_bar, rho_bar, mu_bar, rr_bar);
%         fprintf('   [%10.6f  %10.6f]  [%10.6f  %10.6f]\n', Aa, Ba);

% %         if (a == 10 | a == 30)
% %             if (dumPIdx == -1)
% %                 fprintf('a %d b %d, dumPIdx %d \n',a, b, dumPIdx);
% %             else
% %                 fprintf('a %d b %d, dumPIdx %d \n',a, dummyP.idx(dumPIdx), dumPIdx);
% %             end
% %         end
% %         if (a == 45)
% %             fprintf('a %d b %d, dumPIdx %d \n',a, b, dumPIdx);
% %         end
            
            
    end
%     if (a == 10 | a == 30)
%         fprintf('\n');
%     end
    
%     fprintf('A = %g  %g\n', A);
%     fprintf('B = %g  %g\n', B);
    
    f(:, a) = -A + B + [pb.F;0];
    g(a) = C;
    
  
%     if ((a == 10) | (a == 50))
%         fprintf('value of f %f , %f , %d\n',f(1,a), f(2,a),a);
% %         fprintf('value of g %f , %d\n',g(a),a);
% %         fprintf('basic SPH sum %f, %f  counterV %d  a  %d \n',D(1), D(2), counterV, a);
%     end
    
        
end
% for a = 10:10:100
%     fprintf('value of f %f , %f , %d\n',f(1,a), f(2,a), a);
% end


% Reshape 'f' and 'g' as column vectors.
f = reshape(f, 2 * pb.N, 1);
g = g';

