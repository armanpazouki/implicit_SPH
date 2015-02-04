function dummyP = update_ghosts(edgeP, dummyP)
dummyP.v = edgeP.v(dummyP.idx,:);
dummyP.p = edgeP.p(dummyP.idx);
end
    
