classdef MyFlowSystem < handle
    properties
        pb;
        part;
        edgeP;
        dummyP;
        partNew;
        edgePNew;

        nb_p;
        nb_d;
        %   dr/dt = v; dv/dt = f; g=0 (constraints)
        f; 
        g;
        c;
        %   Jacobian
        Jf; %   Jacobian associated with f [fr fv fp]
        Jg; %   Jacobian associated with g [gr gv gp]
        Jc; %   Jacobian associated with c [gv gp]
        %   jocobian components
        fr;
        fv;
        fp;
        gr;
        gv;
        
        fv_edgeP;
        fp_edgeP;
        gv_edgeP;
        %   Solver Property
        tau;
        beta0;
        %   number of variables
        numPos;
        numVel;
        numPres;
        numVel_edge;
        numPres_edge;
        
        cols_r;
        cols_v;
        cols_p;
        cols_v_edge;
        cols_p_edge;
    end
    methods
        
        function obj = MyFlowSystem()
        end
        
        function InitializeProblem(obj)
            obj.pb = init_problem();
            obj.beta0 = 1;
            obj.tau = obj.pb.dt;
            obj.part = init_particles(obj.pb);  
            [obj.edgeP, obj.dummyP] = set_ghosts(obj.pb, obj.part);
            obj.partNew = obj.part;
            obj.edgePNew = obj.edgeP;
            obj.pb.Ne = obj.edgeP.num;
            obj.pb.NT = obj.pb.N + obj.pb.Ne;
            
            obj.nb_p = cell(obj.pb.NT, 1);
            obj.nb_d = cell(obj.pb.NT, 1);
        end
        
        function dummyP = ApplyBoundaryInternal(obj, edgeP, dummyP)
            dummyP = update_ghosts(edgeP, dummyP);
        end
        
%         function ApplyBoundaryToOldParticles(obj)
%             obj.ghost = obj.ApplyBoundaryInternal(obj.pb, obj.part);
%         end
        
        function ApplyBoundary(obj)            
            obj.dummyP = obj.ApplyBoundaryInternal(obj.edgePNew, obj.dummyP);
        end
        
        function FindNeighbours(obj)
            for i = 1 : obj.pb.NT
                ri = obj.grabR(i, obj.partNew, obj.edgePNew);
                [nb_p1, nb_d1] = find_neighbours(ri, obj.pb, cat(2, obj.partNew.r, obj.edgePNew.r), obj.dummyP);
                obj.nb_p{i} = nb_p1;
                obj.nb_d{i} = nb_d1;
            end
        end
                   
        function CalcRHS(obj)
            [obj.f, obj.g] = rhs(obj.pb, obj.partNew, obj.edgePNew, obj.dummyP);
        end
        
        function CalcJacobian(obj)
            obj.numPos = 2 * obj.pb.N;
            obj.numVel = 2 * obj.pb.N;
            obj.numPres = obj.pb.N;
            obj.numVel_edge = 2 * obj.edgeP.num; %ghost
            obj.numPres_edge = obj.edgeP.num;
            
            obj.cols_r = (1:obj.numPos);
            obj.cols_v = (obj.numPos+1:obj.numPos+obj.numVel);
            obj.cols_p = (obj.numPos+obj.numVel+1:obj.numPos+obj.numVel+obj.numPres);
            domainDataSize = obj.numPos+obj.numVel+obj.numPres;
            obj.cols_v_edge = (domainDataSize + 1 : domainDataSize+obj.numVel_edge);
            obj.cols_p_edge = (domainDataSize + obj.numVel_edge + 1 : domainDataSize + obj.numVel_edge + obj.numPres_edge);
            
            % calc Jacobian
            [obj.Jf, obj.Jg] = jac(obj, obj.pb, obj.partNew, obj.edgePNew);
            % calc sub-Jacobians
            obj.fr = obj.Jf(:,obj.cols_r);
            obj.fv = obj.Jf(:,obj.cols_v);
            obj.fp = obj.Jf(:,obj.cols_p);
            obj.gr = obj.Jg(:,obj.cols_r);
            obj.gv = obj.Jg(:,obj.cols_v);  
            
            obj.fv_edgeP = obj.Jf(:,obj.cols_v_edge);
            obj.fp_edgeP = obj.Jf(:,obj.cols_p_edge);
            obj.gv_edgeP = obj.Jg(:,obj.cols_v_edge);   
        end
        
        function CalcJacobianAndRHSConstraints(obj)
            [obj.Jc, obj.c] = jacAndRHSConstraint(obj, obj.pb, obj.partNew, obj.edgePNew);
        end
        
        function error = Iterate(obj)
            A1 = cat(2, eye(obj.numVel) - (obj.tau * obj.beta0)^2 * obj.fr - (obj.tau * obj.beta0) * obj.fv, ...
                -obj.tau * obj.beta0 * obj.fp,...
                - (obj.tau * obj.beta0) * obj.fv_edgeP, -obj.tau * obj.beta0 * obj.fp_edgeP);
            A2 = cat(2, obj.tau * obj.beta0 * obj.gr + obj.gv, zeros(size(obj.g,1), obj.numPres),...
                obj.gv_edgeP, zeros(size(obj.g,1), obj.numPres_edge));
            Au1 = cat(1, A1, A2);
            A3 = obj.Jc(: , obj.numPos+1:obj.numPos+obj.numVel+obj.numPres+obj.numVel_edge+obj.numPres_edge);
            'yo'
            obj.edgePNew.num


            A = cat(1, Au1, A3);
            size(A);
%             s = svd(A);
%             sizeS = size(s,1);
%             (s([sizeS-10:sizeS]))';
            resR = obj.MakeVertical_2N(obj.partNew.r-obj.tau*obj.beta0*obj.partNew.v-obj.part.r, obj.pb.N);
            b1 = - cat(1, obj.MakeVertical_2N(obj.partNew.v-obj.part.v, obj.pb.N)-obj.tau*obj.beta0*obj.f, obj.g) + ...
                cat(1, -obj.tau * obj.beta0 * obj.fr*resR, obj.gr*resR); 
            b = cat(1, b1, -obj.c);
            %% Add extra pressure constraint
%             myConstraint = zeros(1, obj.numPres + obj.numVel);
%             myConstraint(1, obj.numPres + obj.numVel) = 1;
%             A = cat(1, A, myConstraint);
%             b = cat(1, b, 0);
            %%
            fprintf('size and rank of A %d %d\n', size(A,1), rank(A));
            %% Solve using inverse (assuming full rank)
%             fprintf('rank and size of A %d %d\n', rank(A), size(A,1));
%             res = A\b;
            %% Solve using pinv
            res = pinv(A) * b;
            %% Solve using QR of A'
%             n = size(A,1);
%             [Q, R] = qr(A');
%             % the rank of R should be 2...
%             r = rank(R)
%             % solve for z from R'*z = b
%             % this should really be done with a forward substitution...
%             RR = R(1:r,1:r)';
%             z = zeros(n,1);
%             z(1:r) = RR\b(1:r);
%             % recover solution
%             % this solution has minimum norm (should be same as x2)
%             res = Q * z;
            %% Solve using QR of A
%             n = size(A,1);
%             [Q, R] = qr(A);
%             % the rank of R should be 2...
%             r = rank(R)
%             % modify RHS
%             bb = Q'*b;
%             % solve for x from R*x = bb
%             % this should really be done with a backward substitution
%             % the resulting solution has a minimum number of non-zero elements.
%             res = zeros(n,1);
%             res(1:r) = R(1:r,1:r)\bb(1:r);
            %%
            error = max(abs(res));%/ max(max(max(abs(obj.part.v))) , max(obj.part.p));
            obj.partNew.v = obj.partNew.v + obj.MakeHorizontal_2N(res(1:obj.numVel), obj.pb.N);
            obj.partNew.p = obj.partNew.p + obj.MakeHorizontal_N(res(obj.numVel+1:obj.numVel+obj.numPres), obj.pb.N);
            obj.partNew.r = obj.partNew.r + obj.MakeHorizontal_2N(-resR  + obj.tau * obj.beta0 * res(1:obj.numVel), obj.pb.N);
            
            obj.edgePNew.v = obj.edgePNew.v + ...
                obj.MakeHorizontal_2N(res(obj.numVel+obj.numPres+1 : obj.numVel+obj.numPres+obj.numVel_edge), obj.edgePNew.num);
            obj.edgePNew.p = obj.edgePNew.p + ...
                obj.MakeHorizontal_N(res(obj.numVel+obj.numPres+obj.numVel_edge+1: obj.numVel+obj.numPres+obj.numVel_edge+obj.numPres_edge), obj.edgePNew.num);
        end
        
        function CopyNewToCurrent(obj)
            obj.part = obj.partNew;
            obj.edgeP = obj.edgePNew;
        end
        
        function PeriodicBoundary(obj)
            idx_max = find(obj.partNew.r(1,:) > obj.pb.Lx);
            obj.partNew.r(1,idx_max) = obj.partNew.r(1,idx_max) - obj.pb.Lx;
            idx_min = find(obj.partNew.r(1,:) < 0);
            obj.partNew.r(1,idx_min) = obj.partNew.r(1,idx_min) + obj.pb.Lx;
        end 
        
%         function c = Cols_r(obj, i, n)
%             c = ((i-1)*2 + 1 : i*2);
%         end
% 
%         function c = Cols_v(obj, i, n)
%             numPos = 2 * n;
%             c = (numPos + (i-1)*2 + 1 : numPos + i*2);
%         end
% 
%         function c = Cols_p(obj, i, n)
%             c = 2 * 2 * n + i;
%         end

        function c = Cols_r(obj, i)
            np = obj.pb.N;
            if (i <= np)
                c = ((i-1)*2 + 1 : i*2);
            else
                'something is wrong probably'
                c = [];
            end
        end

        function c = Cols_v(obj, i)
            np = obj.pb.N;
            ne = obj.pb.Ne;
            if (i <= np)
                numSavedComps = 2 * np;
            else
                numSavedComps = 5 * np;
            end
            c = (numSavedComps + (i-1)*2 + 1 : numSavedComps + i*2);
        end
        
        function c = Cols_p(obj, i)
            np = obj.pb.N;
            ne = obj.pb.Ne;
            if (i <= np)
                numSavedComps = 4 * np;
            else
                numSavedComps = 5 * np + 2 * ne; %everything for parts + vel of edgePs
            end
            c = numSavedComps + i;
        end
        
        
        
        
%         function c = Cols_v_edge(obj, i, np, ne)
%             numParPosVelPre = 5 * np;
%             c = (numParPosVelPre + (i-1)*2 + 1 : numParPosVelPre + i*2);
%         end
% 
%         function c = Cols_p_edge(obj, i, np, ne)
%             numParPosVelPre = 5 * np;
%             c = numParPosVelPre + 2 * ne + i;
%         end
        
        % re-shapes the array of (2,N) to (2N,1)
        function out = MakeVertical_2N(obj, r, n)
            if ((size(r, 1) ~= 2) | (size(r, 2) ~= n))
                size(r)
                fprintf('you are calling wrong function: your size is not (2,N)');
            end
            out = zeros(2 * n, 1);
            out(1:2:2 * n) = r(1,:);
            out(2:2:2 * n) = r(2,:);
        end
        
        % re-shapes the array of (1,N) to (N,1)
        function out = MakeVertical_N(obj, r, n)
            if ((size(r, 1) ~= 1) | (size(r, 2) ~= n))
                size(r)
                fprintf('you are calling wrong function: your size is not (1,N)\n');
            end
            out = r';
        end  
        
        % Undo MakeVertical_2N, i.e. re-shape from (2N,1) to (2,N)
        function out = MakeHorizontal_2N(obj, r, n)
            if ((size(r, 1) ~= 2*n) | (size(r, 2) ~= 1))
                size(r, 1)
                fprintf('error in number of components in horizontal vec: not (2N, 1)\n');
            end
            out = zeros(2, n);
            out(1, :) = r(1:2:2*n);
            out(2, :) = r(2:2:2*n);
        end
        
        % Undo MakeVertical_N, i.e. re-shape from (N,1) to (1,N)
        function out = MakeHorizontal_N(obj, r, n)
            if ((size(r, 1) ~= n) | (size(r, 2) ~= 1))
                fprintf('error in number of components in horizontal vec: not (N, 1)');
            end
            out = r';
        end
        
        function [Gr, Gv, Gp] = ghost_influence(obj,bc)
            % The 2 dimensional array 'bc' encodes the BC for some ghost point. We can
            % have the following cases:
            %   bc(1) = 0  or  1
            %   bc(2) = 0  or  2

            if bc(2) == 0
               Gr = [1 0; 0 1];
               Gv = [1 0; 0 1];
            else
               Gr = [1 0; 0 -1]; 
               Gv = [-1 0; 0 -1];
            end

            Gp = 1;
        end
        
        function r = grabR(obj, i, part, edgeP)
            if (i <= obj.pb.N)
                r =part.r(:,i);
            elseif (i <= obj.pb.NT)
                r = edgeP.r(:, i - obj.pb.N);
            else
                'index i out of bound'
            end
        end
        function v = grabV(obj, i, part, edgeP)
            if (i <= obj.pb.N)
                v =part.v(:,i);
            elseif (i <= obj.pb.NT)
                v = edgeP.v(:, i - obj.pb.N);
            else
                'index i out of bound'
            end
        end
        function p = grabP(obj, i, part, edgeP)
            if (i <= obj.pb.N)
                p =part.p(:,i);
            elseif (i <= obj.pb.NT)
                p = edgeP.p(:, i - obj.pb.N);
            else
                'index i out of bound'
            end
        end
        
        function out = IsEdgeP(obj, b) 
            if (b > obj.pb.NT)
                'wrong indexing, IsEdgeP'
                out = 1 > 2;
            end
            out = b > obj.pb.N;
        end
    end %methods
end % classdef
            
            
            
        
        