function Cuts(varargin)
%% debug
% hF = figure(1);
% hAx_Plot = 'New';
% Direction = 'X';

%% Select Figure to be cut, even without the hF input
if eq(nargin,0)
    hF_All = get(0,'Children');
    N_F_All = length(hF_All);

    if eq(N_F_All,0)
        disp('No figures to be cut...')
        return
    end

    AvaliableFigureStr = cell(N_F_All,1);

    for i = 1:N_F_All
        FigName = num2str(hF_All(i).Number);
        AvaliableFigureStr{i} = ['Figure: ',FigName];
    end

    [Choise_hF_ind,~] = listdlg('PromptString','Select the Figure to apply cut:',...
        'SelectionMode','Single',...
        'ListSize',[300,160],...
        'ListString',AvaliableFigureStr);
    hF = hF_All(Choise_hF_ind);
else
    hF = varargin{1};
end
%% Select the axes within the figure and find the X,Y,Z data
hContour = findobj(hF,'type','contour');
          
N_Contour = length(hContour);
switch N_Contour
    case 0
        disp('The selected figure has no contour axis for Cuts...')
        return
    case 1
        Axis_Ind = 1;
    otherwise
        % if more than two avalibale axis list them for slection
        AvaliableAxesStr = cell(N_Contour,1);
        for i = 1:N_Contour
            Title = hContour(i).Parent.Title.String;
            if iscell(Title)
                Title = Title{1};
            end
            AvaliableAxesStr{i} = Title;
        end
        [Axis_Ind,~] = listdlg('PromptString','Select a axis to be cut:',...
                             'SelectionMode','Single',...
                             'ListSize',[300,160],...
                             'ListString',AvaliableAxesStr);
        
end

hContour    = hContour(Axis_Ind);
hAx_Contour = hContour.Parent;

%% Select a figure to place cuts
% get the list of figure
hF_All = get(0,'Children');
N_F_All = length(hF_All);
AvaliableFigureStr = cell(1,1);
AvaliableFigureStr{1} = 'Create New';
for i = 1:N_F_All
    FigName = num2str(hF_All(i).Number);
    AvaliableFigureStr{i+1} = ['Figure: ',FigName];
end

% ask for selection
[Choise_hF_Plot,~] = listdlg('PromptString','Select a figure to place cuts',...
    'SelectionMode','single',...
    'ListSize',[300,160],...
    'ListString',AvaliableFigureStr);

% assign axes
if eq(Choise_hF_Plot,1)
    hF_Cut  = figure;
    hAx_Cut = axes('Parent',hF_Cut);
else
    hF_Cut  = hF_All(Choise_hF_Plot-1);
    hAx_Cut = findobj(hF_Cut,'Type','Axes');
end

%% draw cuts
% add initial imline on the contour axis
XLim   = hAx_Contour.XLim;
X_L    = XLim(2) - XLim(1);
Cut0_X = [XLim(1) + 0.2*X_L, XLim(2) - 0.2*X_L];
Cut0_Y = [1,1].*mean(hAx_Contour.YLim);
hIntL  = imline(hAx_Contour,Cut0_X,Cut0_Y);

% Draw the initial cut on the cut presenting axis
hold(hAx_Cut,'on')
    FakeLineX = linspace(Cut0_X(1),Cut0_X(2));
    FakeLineY = zeros(size(FakeLineX));
    hCutLine = plot(hAx_Cut,FakeLineX,FakeLineY,'LineWidth',2);
    CutCallback([Cut0_X;Cut0_Y]',hContour,hCutLine)
hold(hAx_Cut,'off')

hAx_Cut.XGrid = 'on';
hAx_Cut.YGrid = 'on';
hAx_Cut.XLim  = hAx_Contour.XLim;
ZMax = max(hContour.ZData(:));
hAx_Cut.YLim   = [-ZMax,ZMax].*1.1;
hAx_Cut.XLabel.String = hAx_Contour.XLabel.String;
hAx_Cut.YLabel.String = 'Z intensity';
hAx_Cut.FontSize = 16;

OriginalTitle = hAx_Contour.Title.String;
if iscell(OriginalTitle)
    OriginalTitle = OriginalTitle{1};
end
hAx_Cut.Title.String = ['X-Z cut of ',OriginalTitle];
setColor(hIntL,hCutLine.Color)

addNewPositionCallback( hIntL,@(pos) CutCallback(pos,hContour,hCutLine) );

function CutCallback(pos,hContour,hCutLine)
%% Find the X,Y,Z data
X = hContour.XData;
Y = hContour.YData;
Z = hContour.ZData;

%% Generate Cut value Vq
% pos = [X1,Y1;X2,Y2]
Xq = linspace(pos(1),pos(2));
Yq = linspace(pos(3),pos(4));

Vq = interp2(X,Y,Z,Xq,Yq);

%% Update the cut line on hCut
hCutLine.XData = Xq;
hCutLine.YData = Vq;
