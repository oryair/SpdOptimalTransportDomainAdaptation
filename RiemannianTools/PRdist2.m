function M = PRdist2(PP1, PP2)
    N1 = length(PP1);
    N2 = length(PP2);
    M  = zeros(N1, N2);
    
    for ii = 1 : N1
        for jj = 1 : N2
            M(ii,jj) = RiemannianDist(PP1{ii}, PP2{jj});
        end
    end
end
