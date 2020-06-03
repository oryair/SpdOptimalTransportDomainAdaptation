function mLamP = SinkhornRegOptimalTransport(Covs1, Covs2, vClass1, lambda)

    KDE = true;
    N1  = length(Covs1);
    N2  = length(Covs2);
    
    if KDE == true
        mWx  = CalcDist(Covs1);
        epsX = .1 * median(mWx(:));
        mKx  = exp(-mWx.^2 / (2 * epsX^2));
        vP1  = sum(mKx, 2);
        vP1  = vP1 / sum(vP1);
        
        mWy  = CalcDist(Covs2);
        epsY = .1 * median(mWy(:));
        mKy  = exp(-mWy.^2 / (2 * epsY^2));
        vP2  = sum(mKy, 2);
        vP2  = vP2 / sum(vP2);
    else
        vP1  = ones(N1, 1) / N1;
        vP2  = ones(N2, 1) / N2;
    end
    
    M  = CalcDist2(Covs1, Covs2).^2;
    
    if nargin < 4
        sig    = 0.05 * median(M(:));
        lambda = 1 / (2 *  sig^2);
    end
       
    %--
    G  = zeros(size(M));
    p  = 1/2;
    vC = unique(vClass1);
    Nc = length(vC);

    
    numIter = 10;
    for kk = 1 : numIter
        M = M + G;
        
        %--
        K    = exp(-lambda * M);
        L    = length(vP1);
        u    = ones(L, 1) / L;
        K2   = K ./ vP1;
        
        for ii = 1 : 1000
            u  = 1 ./ (K2 * (vP2 ./ (K' * u)));
        end
        v    = vP2 ./ (K' * u);
        mLamP = (K .* u) .* v';
        
        %--
        for jj = 1 : N2
            for cc = 1 : Nc
                vIdx = vClass1 == vC(cc);
                G(vIdx,jj) = p * (norm(mLamP(vIdx,jj), 1) + 1e-6)^(p-1);
            end
        end
    end
end

%%
function mW = CalcDist(Covs)
    N  = length(Covs);
    mW = zeros(N, N);
    
    for ii = 1 : N
        for jj = ii + 1 : N
            mW(ii,jj) = RiemannianDist(Covs{ii}, Covs{jj});
        end
    end
    mW = mW + mW';
end

%%
function M = CalcDist2(Covs1, Covs2)
    N1 = length(Covs1);
    N2 = length(Covs2);
    M  = zeros(N1, N2);
    
    for ii = 1 : N1
        for jj = 1 : N2
            M(ii,jj) = RiemannianDist(Covs1{ii}, Covs2{jj});
        end
    end
end