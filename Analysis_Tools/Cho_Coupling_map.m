% Plot_Cho_Coupling_Matrix
Beta = hMain.FTIR.H.Beta;
hF = figure;
pcolor(Beta);
hAx = findobj(hF,'Type','axes');

% colormap
colorbar
hAx.CLim = [-12,12];
Cho_cmap = load('Cho_cmap');
colormap(hF,Cho_cmap.Cho_cmap)
