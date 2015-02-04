function edgeP = update_edgeP(obj, part, edgeP, dummyP, pb)
%% set velocity
for a = 1 : pb.Ne
    bcType = edgeP.bc(:, a);
    if (bcType(1) == -2) %left
       edgeP.v(:, a) = [pb.uL;pb.vL];
    elseif (bcType(1) == 2) %right
       edgeP.v(:, a) = [pb.uR;pb.vR];
    elseif (bcType(2) == 2)
       edgeP.v(:, a) = [pb.uT;pb.vT];
    else %(bcType(2) == -2)
       edgeP.v(:, a) = [pb.uB;pb.vB];    
    end
end
%% calc pressure
for ia = 1 : pb.Ne
    a = ia + pb.N;  % refer to edge index in the after collision index array
    nb_p = obj.nb_p{a};
    nb_d = obj.nb_d{a};
    ra = obj.grabR(a, part, edgeP);
    pa_sph = 0;
%     pa_sph =  pb.m / pb.rho * obj.grabP(a, part, edgeP) * kernel(zeros(size(ra)), pb.h, 0);
    for ib = 1 : length(nb_p) + length(nb_d)
        if (ib <= length(nb_p)) % part or edgeP
            b = nb_p(ib);
%             if (b == a)
%                 continue;
%             end
            r = obj.grabR(a, part, edgeP) - obj.grabR(b, part, edgeP);
        else                    % dummyP
            dumPIdx = nb_d(ib - length(nb_p)); % index of dummy particle in dummyP
            b = dummyP.idx(dumPIdx); % associated edgeP particle  
            b = b + pb.N;            % to comply with the grabV and grabP functions
        	r = obj.grabR(a, part, edgeP) - dummyP.r(:,dumPIdx);
        end
        pa_sph = pa_sph + pb.m / pb.rho * obj.grabP(b, part, edgeP) * kernel(r, pb.h, 0);
    end
    edgeP.p(:, ia) = pa_sph;
%%
end
    
