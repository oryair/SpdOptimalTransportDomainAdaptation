function mGamma = SinkhornOptimalTransport(mC, vP1, vP2, lambda)
  
    N1 = length(vP1);
%     N2 = length(vP2);

    if nargin < 4
        sig    = 0.05 * median(mC(:));
        lambda = 1 / (2 *  sig^2);
    end
    
    K  = exp(-lambda * mC);
    u  = ones(N1, 1) / N1;
    K2 = K ./ vP1;
    
    for ii = 1 : 1000
        u  = 1 ./ (K2 * (vP2 ./ (K' * u)));
    end
    
    v      = vP2 ./ (K' * u);
    mGamma = (K .* u) .* v';
end
