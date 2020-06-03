close all
clear

addpath('../RiemannianTools/');
addpath(genpath('../OptimalTransportTools/'));

rng(1);

%%
N = 50;
D = 2;

Covs1{N} = [];
Covs2{N} = [];
        
% T = [0.9 .5;
%     .5    1.1];

% T = [.5   -.25;
%      -.25  1];

T = RandP(D);
 
Nt  = 201;
vTh = linspace(0, pi, Nt);
vR  = nan(Nt, 1);

for ii = 1 : N
    Covs1{ii} = RandP(D);
end

vP1 = ones(N, 1) / N;
vP2 = ones(N, 1) / N;

figure;
for tt = 1 : Nt

    th = vTh(tt);
    U  = [cos(th), sin(th);
         -sin(th), cos(th)];
    
    for ii = 1 : N
        Covs2{ii} = T * U * Covs1{ii} * U' * T';
    end

    mC        = PRdist2(Covs1, Covs2).^2;
    mPlan     = SolveOT(mC, vP1, vP2);
    [vIdx, ~] = find(mPlan);
    
    vD = nan(N, 1);
    for ii = 1 : N
        vD(ii) = RiemannianDist(Covs2{vIdx(ii)}, Covs2{ii})^2;
    end
    vR(tt) = sqrt(mean(vD));
    
    plot(vTh, vR, 'b'); drawnow;
    
%     figure; imagesc(mPlan); colorbar;
%     keyboard
end
close;

%%
figure; hold on; grid on; set(gca, 'FontSize', 14);
plot(vTh, vR, 'b', 'LineWidth', 3);
xlabel('$\theta$', 'Interpreter', 'latex');
% title('Correct matching percentage', 'Interpreter', 'latex');
h = legend('$\sqrt{\frac{1}{N}\sum_{i=1}^{N}d_{R}^{2}\left(t_{\theta} \left(P_{i}\right),S_{\theta}P_{i}S_{\theta}^{T}\right)}$');
set(h, 'Interpreter', 'latex');
axis tight

%%
vTh  = [0, pi/20, pi/10, pi/6];
cStr = {'$0$', '$\frac{\pi}{20}$', '$\frac{\pi}{10}$', '$\frac{\pi}{6}$'};

for tt = 1 : length(vTh)

    th = vTh(tt);
    U  = [cos(th), sin(th);
         -sin(th), cos(th)];
    
    for ii = 1 : N
        Covs2{ii} = T * U * Covs1{ii} * U' * T';
    end

    mC    = PRdist2(Covs1, Covs2).^2;
    mPlan = SolveOT(mC, vP1, vP2);
    figure; imagesc(mPlan); % colorbar;
    axis equal; axis tight; set(gca, 'FontSize', 16);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function P = RandP(D)
    P = randn(D);
    P = P * P' + .5 * eye(D);
end

