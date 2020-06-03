function mGamma = SolveOT(mC, vP1, vP2)

if nargin == 1
    [N1, N2] = size(mC);
    vP1      = ones(N1, 1) / N1;
    vP2      = ones(N2, 1) / N2;
end

flat = @(x) x(:);

Cols = @(n0, n1) sparse( flat(repmat(1:n1, [n0 1])), ...
                         flat(reshape(1:n0*n1,n0,n1) ), ...
                         ones(n0*n1,1) );

Rows = @(n0, n1) sparse( flat(repmat(1:n0, [n1 1])), ...
                         flat(reshape(1:n0*n1,n0,n1)' ), ...
                         ones(n0*n1,1) );

Sigma = @(n0, n1) [Rows(n0,n1);
                   Cols(n0,n1)];
              
maxit   = 5e4;
tol     = 1e-9;
otransp = @(C,p0,p1) reshape(perform_linprog(                   ...
                                  Sigma(length(p0),length(p1)), ...
                                  [p0(:); p1(:)],               ...
                                  C(:), 0, maxit, tol),         ...
                                  [length(p0) length(p1)] );
                  

mGamma = otransp(mC, vP1, vP2);

end