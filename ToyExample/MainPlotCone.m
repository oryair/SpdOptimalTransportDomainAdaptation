close all
clear

addpath('../RiemannianTools/');
addpath(genpath('../OptimalTransportTools/'));

%%
N = 5;
D = 2;

Covs1{N} = [];
Covs2{N} = [];
        
T = [.5   -.25;
     -.25  1];
 
% T = [0.9 -.5;
%      -.5    1.1];

Covs1{1} = [6  2;
            2  3];
Covs1{2} = [4    1.5;
            1.5  6];
Covs1{3} = [3.5  2;
            2    5];
Covs1{4} = [6    1;
            1    5];
Covs1{5} = [7  2;
            2  3.5];

vTh(1) = 0;      
vTh(2) = pi / 2; 

legendStrS{1} = '$S_{0}P_{i}S_{0}^{T}$';
legendStrS{2} = '$S_{\frac{\pi}{2}}P_{i}S_{\frac{\pi}{2}}^{T}$';

legendStrT{1} = '$t_0\left(P_{i}\right)$';
legendStrT{2} = '$t_{\frac{\pi}{2}}\left(P_{i}\right)$';

mColor = lines(N);

for tt = 1 : 2
    
    th  = vTh(tt);
    U   = [cos(th), sin(th);
          -sin(th), cos(th)];
    
    S = T * U;
    for ii = 1 : N
        Covs2{ii} = S * Covs1{ii} * S';
    end
    
    figure; hold on; set(gca, 'FontSize', 16);
    PlotCone(12);
    view([-110, 9]);
    axis tight
    
    for ii = 1 : N
        PlotPair(Covs1{ii}, Covs2{ii}, mColor(ii,:), legendStrS{tt});
    end
    
    mC        = PRdist2(Covs1, Covs2).^2;
    mPlan     = SolveOT(mC);
    [vIdx, ~] = find(mPlan);
    Covs1OT   = Covs2(vIdx);
    
    figure; hold on; set(gca, 'FontSize', 16);
    PlotCone(12);
    view([-110, 9]);
    axis tight
    
    for ii = 1 : N
        PlotPair(Covs1{ii}, Covs1OT{ii}, 'm', legendStrT{tt});
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPair(P1, P2, c, str)
    h1 = plot3(P1(1), P1(2), P1(4), 'b.', 'MarkerSize', 40);
    h2 = plot3(P2(1), P2(2), P2(4), 'r.', 'MarkerSize', 40);
    
    t = linspace(0, 1, 100);
    for ii = 1 : 100
        tt = t(ii);
        P = P1^(1/2) * (P1^(-1/2) * P2 * P1^(-1/2))^tt * P1^(1/2);
        vX(ii) = P(1);
        vY(ii) = P(2);
        vZ(ii) = P(4);
    end
    plot3(vX, vY, vZ, 'Color', c, 'LineWidth', 3);
  
    h = legend([h1, h2], '$P_i$', str, 'Location', 'northwest'); set(h, 'Interpreter', 'Latex');
end

%%
function PlotP(P, c)
    plot3(P(1), P(2), P(4), [c, '.'], 'MarkerSize', 20);
end