function [Covs1OT] = ApplyPlan(Covs1, Covs2, mPlan, K)

    N1 = length(Covs1);

    if nargin < 4
        K = 10;
    end
    K = min(K, N1);
    
    [mW, mIdx] = sort(mPlan', 'descend');
    
    mW         = mW(1:K,:);
    mW         = mW ./ sum(mW);
    mIdx       = mIdx(1:K,:);
    
    N1      = length(Covs1);
    Covs1OT = Covs1;
    for ii = 1 : N1
        Covs1OT{ii} = WeightedRiemannianMean(cat(3, Covs2{mIdx(:,ii)}), mW(:,ii));
    end
end
