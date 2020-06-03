function PlotCone(R)

    N = 50;
    x = linspace(0,  R, N);
    z = linspace(0,  R, N);
    y = linspace(-R, R, N);
    
    [XX, YY, ZZ] = meshgrid(x, y, z);
    mX = [XX(:), YY(:), ZZ(:)]';
    
    vIdx1 = mX(2,:).^2 < mX(1,:) .* mX(3,:);
    vIdx2 = mX(1,:) + mX(3,:) < R;
    vIdx  = vIdx1 & vIdx2;

    vIdx2 = find(vIdx);
    L     = length(vIdx2);
    vC    = nan(L, 1);
    for ll = 1 : L
        M      = reshape(mX([1, 2, 2, 3],vIdx2(ll)), 2, 2);
        vC(ll) = min(eig(M));
    end
    
    
    K = convhull(mX(1,vIdx), mX(2,vIdx), mX(3,vIdx));
    trisurf(K, mX(1,vIdx), mX(2,vIdx), mX(3,vIdx), 'FaceColor', [0.1, 0.1, 0.1], 'FaceAlpha', 0.05, 'LineStyle', 'None');
%     trisurf(K, mX(1,vIdx), mX(2,vIdx), mX(3,vIdx), vC, 'FaceAlpha', .5, 'LineStyle', 'None');
    
%     h = scatter3(mX(1,vIdx), mX(2,vIdx), mX(3,vIdx), 25, [0.1, 0.1, 0.1], 'Fill');
%     h = scatter3(mX(1,vIdx), mX(2,vIdx), mX(3,vIdx), 250, vC, 'Fill');
%     alpha(h, .05);

    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    zlabel('$z$', 'Interpreter', 'latex');
    axis equal;
    grid on;
    
end