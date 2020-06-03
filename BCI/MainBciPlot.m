close all;
clear;

addpath('../RiemannianTools');
addpath(genpath('../OptimalTransportTools/'));

%% Load Data:
S1 = 8;
D1 = 1;

S2 = 8;
D2 = 2;

disp('Load Data');
[Events1, vClass1] = GetEvents(S1, D1);
Covs1              = CalcCovs(Events1);

[Events2, vClass2] = GetEvents(S2, D2);
Covs2              = CalcCovs(Events2);

vS     = [1 * ones(1, length(vClass1)), 2 * ones(1, length(vClass2))];
vClass = [vClass1, vClass2];

%% No OT
Covs   = [Covs1, Covs2];
mX     = CovsToVecs(cat(3, Covs{:}));

%%
mTSNE = tsne(mX')';
figure; PlotData(mTSNE, vClass, vS);
subplot(1,2,1); title('Before applying OT', 'Interpreter', 'Latex');

%% Applying OT
N1      = length(Events1);
N2      = length(Events2);
vP1     = ones(N1, 1) / N1;
vP2     = ones(N1, 1) / N2;
mC      = PRdist2(Covs1, Covs2).^2;
mPlan   = SinkhornRegOptimalTransport(Covs1, Covs2, vClass1); %-- Supervised
% mPlan   = SinkhornOptimalTransport(mC, vP1, vP2); %-- Unsupervised
Covs1OT = ApplyPlan(Covs1, Covs2, mPlan, 20);

%% 
Covs   = [Covs1OT, Covs2];
mX     = CovsToVecs(cat(3, Covs{:}));

%%
mX1 = mX(:, vS == 1);
mX2 = mX(:, vS == 2);
linaerSvmTemplate = templateSVM('Standardize', false);
mdlLinearSVM      = fitcecoc(mX1', vClass1, 'Learners', linaerSvmTemplate);
res               = mean( mdlLinearSVM.predict(mX2') == vClass2' )

%%
mTSNE = tsne(mX')';
figure; PlotData(mTSNE, vClass, vS);
subplot(1,2,1); title('After applying OT', 'Interpreter', 'Latex');
subplot(1,2,2); title(['Accuracy - $', num2str(100*res), '\%$'], 'Interpreter', 'Latex');


%%
function Covs = CalcCovs(Events)
    for ii = 1 : length(Events)
        mX       = Events{ii}';
        Covs{ii} = cov(mX');
    end
end

%%
function PlotData(mX, vClass, vS)

    vMarker = 'od';
    vColorS = 'br';
    vColorC = 'gmyk';
    
    vUniqueS  = unique(vS);
        
    %--
    subplot(1,2,2);
    for cc = 1 : 4
        vIdxC = vClass == cc;
        for ss = 1 : 2
            marker = vMarker(ss);
            vIdxS  = vS == vUniqueS(ss);
            color  = vColorS(ss);
            vIdx   = vIdxS & vIdxC;
            scatter(mX(1,vIdx), mX(2,vIdx), 50, color, marker, 'Fill', 'MarkerEdgeColor', 'k'); hold on;
        end
    end
    
    h = legend({['Subject - ', num2str(vUniqueS(1))];
                ['Subject - ', num2str(vUniqueS(2))]}, ...
                'FontSize', 12, 'Location', 'Best'); set(h, 'Color', 'None');
    axis tight;
    
    %--
    subplot(1,2,1);
    for ss = 1 : 2
        marker = vMarker(ss);
        vIdxS  = vS == vUniqueS(ss);
        for cc = 1 : 4
            color  = vColorC(cc);
            vIdxC = vClass == cc;
            vIdx  = vIdxS & vIdxC;
            scatter(mX(1,vIdx), mX(2,vIdx), 50, color, marker, 'Fill', 'MarkerEdgeColor', 'k'); hold on;
        end
    end
    
    h = legend({'Left Hand', 'Right Hand', 'Foot', 'Tongue'}, 'FontSize', 12, 'Location', 'Best'); set(h, 'Color', 'None');
    axis tight;
end
