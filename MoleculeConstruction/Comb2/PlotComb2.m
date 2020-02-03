function hF = PlotComb2(hAx,Structure,~)
%% Prep
if ~ishandle(hAx)
    hF = figure; 
    hAx = axes('Parent',hF);
else
    hF = hAx.Parent;
end

% run through all sub-structures
N_S = length(Structure.Children);
for i = 1:N_S
    S = Structure.Children(i);    
    S.SD_Draw(hAx);
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
