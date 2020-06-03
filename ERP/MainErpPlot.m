close all
clear

addpath('../RiemannianTools/');
addpath(genpath('../OptimalTransportTools/'));

dirPath = 'PATH/TO/DATA/';

%%
vSubjectIdx = [15, 16, 18];
Ns          = length(vSubjectIdx);
L           = 480;
Covs{L,Ns}  = [];

for ii = 1 : Ns
    fileName   = [dirPath, 'Subject', num2str(vSubjectIdx(ii)), 'Session1.mat'];
    load(fileName);
    Data       = mX;
    Labels{ii} = vY(1:L);
    
    vIdx      = vY == 1;
    mE        = squeeze( mean(Data(vIdx,1:16,:), 1) );
%     mE        = CalcERP(Data);
    
    for ll = 1 : L
        mXi         = squeeze( Data(ll,1:16,:) );
        mXX         = [mE;
                       mXi];
        Covs{ll,ii} = cov(mXX');
    end
end

%%
vS = repmat(1:Ns, L, 1);
vS = vS(:);
vC = cat(2, Labels{:});

%% No OT
CovsA = cat(3, Covs{:});
mX    = CovsToVecs(CovsA);

%%
mTSNE = tsne(mX')';
figure; PlotData(mTSNE, vC, vS');
set(gcf, 'Position', [376   763   864   335]);
subplot(1,2,1); title('Before applying OT', 'Interpreter', 'Latex');

%%
for ii = 2 : Ns
    mPlan      = SinkhornRegOptimalTransport(Covs(:,ii), Covs(:,1), Labels{ii}, 0.1);
    Covs(:,ii) = ApplyPlan(Covs(:,ii), Covs(:,1), mPlan, 10);
end

%%
CovsB = cat(3, Covs{:});
mX    = CovsToVecs(CovsB);

%%
mTSNE = tsne(mX')';

%%
figure; PlotData(mTSNE, vC, vS');
set(gcf, 'Position', [376   763   864   335]);
subplot(1,2,1); title('After applying OT', 'Interpreter', 'Latex');

%%
function PlotData(mX, vClass, vS)

    vMarker = 'ods';
    vColorS = 'brg';
    vColorC = 'mk';
    
    vUniqueS  = unique(vS);
        
    %--
    subplot(1,2,2);
    for cc = 1 : 2
        vIdxC = vClass == cc - 1;
        for ss = 1 : 3
            marker = vMarker(ss);
            vIdxS  = vS == vUniqueS(ss);
            color  = vColorS(ss);
            vIdx   = vIdxS & vIdxC;
            scatter(mX(1,vIdx), mX(2,vIdx), 50, color, marker, 'Fill', 'MarkerEdgeColor', 'k'); hold on;
        end
    end
    
%     h = legend({['Subject - ', num2str(vUniqueS(1))];
%                 ['Subject - ', num2str(vUniqueS(2))]}, ...
%                 'FontSize', 12, 'Location', 'Best'); set(h, 'Color', 'None');
    h = legend({'Subject 3'; 'Subject 4'; 'Subject 5'}, ...
                'FontSize', 12, 'Location', 'Best'); set(h, 'Color', 'None');
    axis tight;
    
    %--
    subplot(1,2,1);
    for ss = 1 : 3
        marker = vMarker(ss);
        vIdxS  = vS == vUniqueS(ss);
        for cc = 1 : 2
            color  = vColorC(cc);
            vIdxC = vClass == cc - 1;
            vIdx  = vIdxS & vIdxC;
            scatter(mX(1,vIdx), mX(2,vIdx), 50, color, marker, 'Fill', 'MarkerEdgeColor', 'k'); hold on;
        end
    end
    
    h = legend({'Non-target', 'Target'}, 'FontSize', 12, 'Location', 'Best'); set(h, 'Color', 'None');
    axis tight;
end

%%
function mE = CalcERP(Data)
    load ErpDetect.mat
    
    x = permute(Data(:,1:16,:), [2, 3, 1]);
    for ii = 1 : 480
        C         = cov(x(:,:,ii)');
        x(:,:,ii) = sqrtm(inv(C)) * x(:,:,ii);
    end
    xTest = reshape(x,size(x,1)*size(x,2),size(x,3));
    
    p = classify(ErpDetect, xTest);
    
%     vP = find(sum(p) > 80);
    [~, vP] = sort(sum(p), 'descend');
    vP      = vP(1:20);
    mE = squeeze( mean(Data(vP,1:16,:), 1) );
end