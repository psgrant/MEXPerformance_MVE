clf
[Properties] = ComputeMaterialProperties(X,Theta*Tbp,NodeRho0,Constants);

MAX = max(EWPMesh.NodePos);
MIN = min(EWPMesh.NodePos);
AXLIMS = [MIN(1) MAX(1) MIN(2) MAX(2) MIN(3) MAX(3)];
% Water colormap on subplot 1
ax1 = subplot(2,4,[1 2]);
map = watercolmap();
colormap(ax1,map)
colormap(map)
title('Moisture Content')
[NodeCol, CMin, CMax] = ReturnColsFromData((X), map);
% All the triangular faces
p0 = patch('Faces', PlottingTriangles, 'Vertices', NodePos, 'FaceVertexCData', NodeCol, ...
    'FaceColor', 'interp', 'EdgeColor', 'black', 'FaceAlpha', FAlpha, 'EdgeAlpha', EAlpha);
% All the rectangular faces
p1 = patch('Faces', PlottingQuads, 'Vertices', NodePos, 'FaceVertexCData', NodeCol, ...
    'FaceColor', 'interp', 'EdgeColor', 'black', 'FaceAlpha', FAlpha, 'EdgeAlpha', EAlpha);
clim([CMin, CMax])
PlottingGubbins(AXLIMS)
colorbar()

% Temperature colormap on subplot 2
ax2 = subplot(2,4,[3 4]);
map = watercolmap();
colormap(ax2,map)

[NodeCol, CMin, CMax] = ReturnColsFromData(Properties.Sw, map,0,1);
% All the triangular faces
p0 = patch('Faces', PlottingTriangles, 'Vertices', NodePos, 'FaceVertexCData', NodeCol, ...
    'FaceColor', 'interp', 'EdgeColor', 'black', 'FaceAlpha', FAlpha, 'EdgeAlpha', EAlpha);
% All the rectangular faces
p1 = patch('Faces', PlottingQuads, 'Vertices', NodePos, 'FaceVertexCData', NodeCol, ...
    'FaceColor', 'interp', 'EdgeColor', 'black', 'FaceAlpha', FAlpha, 'EdgeAlpha', EAlpha);
clim([CMin, CMax])
PlottingGubbins(AXLIMS)
colorbar()
title('Saturation')


ax3 = subplot(2,4,[6 7]);
map = TempColMap();
colormap(ax3,map)

[NodeCol, CMin, CMax] = ReturnColsFromData(Theta*100, map,20,30);
% All the triangular faces
p0 = patch('Faces', PlottingTriangles, 'Vertices', NodePos, 'FaceVertexCData', NodeCol, ...
    'FaceColor', 'interp', 'EdgeColor', 'black', 'FaceAlpha', FAlpha, 'EdgeAlpha', EAlpha);
% All the rectangular faces
p1 = patch('Faces', PlottingQuads, 'Vertices', NodePos, 'FaceVertexCData', NodeCol, ...
    'FaceColor', 'interp', 'EdgeColor', 'black', 'FaceAlpha', FAlpha, 'EdgeAlpha', EAlpha);
clim([CMin, CMax])
PlottingGubbins(AXLIMS)
colorbar()
title('Temperature')

%%
if n > 2 % Avg Temp and Moisture plot
    subplot(2,4,5)
    yyaxis left
    plot(TimeArray(1:n)/3600,PhiAvgRho0(1:n),'color',[0.4 0.3 1],'linewidth',2)
    yyaxis right
    plot(TimeArray(1:n)/3600,PhiAvgRho0T(1:n),'color',[1 0.3 0.4],'linewidth',2)
    yyaxis left
    title('X and T')
    xlabel('Time [h]')
    
    ylabel('Moisture Content')
    yyaxis right
    % ylim([20 40])
    ylabel('Temperature')
    grid on
    grid minor
    legend('Moisture','Temperature','location','northwest')
    pbaspect([1.1 1 1])

    subplot(2,4,8) % Kyrlov and Dt plot

    yyaxis right
    p2 = plot(TimeArray(1:n)/3600,TnArray(1:n),'color',[1 0.3 0.4],'linewidth',2);
    yyaxis left
    p1 = plot(TimeArray(1:n)/3600,SigmaSave(1:n),'color',[0.4 0.3 1],'linewidth',2);
    
    yyaxis left
    title('\sigma_w and T_n')
    xlabel('Time [h]')
    % ylim([0,MaxKrylov+2])
    % ylabel('Krylov Dimension [m]')
    ylabel('\sigma_w')
    yyaxis right
    % ylim([20 40])
    ylabel('T_n [s]')
    grid on
    % grid minor
    legend('\sigma','T_n','location','northwest')
    set(gca, 'SortMethod', 'depth')
    p1.ZData = ones(size(p1.XData));
    p2.ZData = zeros(size(p2.XData));
    pbaspect([1.1 1 1])
end



function [] = PlottingGubbins(AXLIMS)
view([-35 20]);
pbaspect(AXLIMS(2:2:6)./max(AXLIMS))
axis(AXLIMS)
end