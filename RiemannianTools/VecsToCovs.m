function Covs = VecsToCovs(mX, M)
    
    MP      = sqrtm(M);
    [D2, N] = size(mX);
    D       = round( (sqrt(8*D2 + 1) - 1) / 2 );
    
    Covs{N} = [];
    mIdx    = triu(true(D, D));
    
    mW = sqrt(2) * ones(D) + (2 - sqrt(2)) * eye(D);
    for ii = 1 : N
        Si       = zeros(D, D);        
        Si(mIdx) = mX(:,ii);
        Si       = (Si + Si') ./ mW;
        Covs{ii} = MP * expm(Si) * MP;
    end
   
end