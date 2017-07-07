function hF = PlotComb2(GUI_data)
%% Prep
Structure = GUI_data.Structure;

hF = figure;
hAx = axes('Parent',hF);

% run through all sub-structures
N_S = length(Structure.Children);
for i = 1:N_S
    S = Structure.Children(i);
    hPlotFunc = S.hPlotFunc;
    GUI_Data = feval(S.hParseGUIFunc,S.hGUIs);

    feval(hPlotFunc,hAx,S,GUI_Data);
end

%% figure options
axis(hAx,'image');
rotate3d(hAx,'on')
grid(hAx,'on')
box(hAx,'on')
view(hAx,[30,15])
hAx.XLabel.String = 'X';
hAx.YLabel.String = 'Y';
hAx.ZLabel.String = 'Z';