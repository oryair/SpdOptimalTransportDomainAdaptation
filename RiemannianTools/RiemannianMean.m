function M = RiemannianMean(tC)

Symm = @(M) (M + M') / 2;

Np = size(tC, 3);

if Np == 1
    M = tC;
    return
end

M  = mean(tC, 3);
vW = ones(Np, 1) / Np;

for ii = 1 : 30
    A = sqrtm(M);
    B = inv(A);
        
    S = zeros(size(M));
    for jj = 1 : Np
        C   = tC(:,:,jj);
        BCB = Symm(B * C * B);
        S   = S + vW(jj) * logm(BCB);
    end
    
    M = Symm(A * expm(S) * A); 
    
    eps = norm(S, 'fro');
    if (eps < 1e-6)
        break;
    end
end

end