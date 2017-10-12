function Cuts(hF,Direction)
%% debug
% hF = figure(1);
% hAx_Plot = 'New';
% Direction = 'X';

%% Select Figure to be cut
if ~ishandle(hF)
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
end

%% Select the axes within the figure and find the X,Y,Z data
h2D = findobj(hF   ,'type','surf',...
              '-or','type','contour',...
              '-or','type','mesh');
Choise = 1;  
N_Axes = length(h2D);
if eq(N_Axes,0)
    disp('Input figure is not valid for Cuts...')
    return
end

if gt(N_Axes,1)
    AvaliableAxesStr = cell(N_Axes,1);
    for i = 1:N_Axes
        Title = h2D(i).Parent.Title.String;
        if iscell(Title)
            Title = Title{1};
        end
        AvaliableAxesStr{i} = Title;
    end
    [Choise,~] = listdlg('PromptString','Select the axes to be cut:',...
        'SelectionMode','Single',...
        'ListSize',[300,160],...
        'ListString',AvaliableAxesStr);
end

h2D = h2D(Choise);
hAx_Parent = h2D.Parent;
X = h2D.XData;
Y = h2D.YData;
Z = h2D.ZData;

%% rename variables
XLim   = hAx_Parent.XLim;
YLim   = hAx_Parent.YLim;
X_Cut  = mean(XLim);
Y_Cut  = mean(YLim);
XLabel = hAx_Parent.XLabel.String;
YLabel = hAx_Parent.YLabel.String;

[~,Slice_X_Ind] = min(abs(X - X_Cut));
[~,Slice_Y_Ind] = min(abs(Y - Y_Cut));

% Env_H =  abs(Z(Slice_Y_Ind,:));
Osc_H = real(Z(Slice_Y_Ind,:));
% Env_V =  abs(Z(:,Slice_X_Ind));
Osc_V = real(Z(:,Slice_X_Ind));

ZMax = max(abs(Z(:)));

%% creat axes for accomodating imline if needed
% if the plot method is other than contour, then I need to create another
% axes on top of the original axis to accomodate imlines

Type_2D = h2D.Parent.Tag;
switch Type_2D
    case 'contour'
        hAx_imline = hAx_Parent;
    otherwise
        hOld_imline = findobj(hF,'Tag','imline');
        assignin('base','hOld_imline',hOld_imline)
        if ishandle(hOld_imline)
            % if there is an axes deficated for imline, then use the old one
            hAx_imline = findobj(hOld_imline,'Type','axes');
        else
            % if not, create a new one.
            hAx_imline = axes('Parent',hF);
            hAx_imline.Tag = 'imline';
            hAx_imline.Position = hAx_Parent.Position;
            hAx_imline.Color = 'non';
            hAx_imline.XLim = XLim;
            hAx_imline.YLim = YLim;
            hAx_imline.XTick = [];
            hAx_imline.YTick = [];
            hAx_imline.Box = 'off';
        end
end

%% Decide cutting direction
switch Direction
    case 'X'
        CursorColor = [255,105,180]./256;
        Osc_X    = X;
        Osc_Y    = Osc_H;
        Osc_C    = [1,0,0];
        
        dX = abs(diff(XLim));
        IntL_X   = [XLim(1)-dX,XLim(2)+dX];
        IntL_Y   = [Y(Slice_Y_Ind),Y(Slice_Y_Ind)];    
    case 'Y'        
        CursorColor = [0,1,0];
        
        Osc_X    = Y;
        Osc_Y    = Osc_V;
        Osc_C    = [0,0,1];
        
        dY = abs(diff(YLim));
        IntL_X   = [X(Slice_X_Ind),X(Slice_X_Ind)];
        IntL_Y   = [YLim(1)-dY,YLim(2)+dY];
end

%% Make figure
% Prep axes if it wasn't assigned
% hAx_Plot = 'New'; % Assume to create a new figure for cuts

% get the list of figure
hF_All = get(0,'Children');
N_F_All = length(hF_All);
AvaliableFigureStr = cell(N_Axes+1,1);
AvaliableFigureStr{1} = 'Create New';
for i = 1:N_F_All
    FigName = num2str(hF_All(i).Number);
    AvaliableFigureStr{i+1} = ['Figure: ',FigName];
end

% ask for selection
[Choise_hF_Plot,~] = listdlg('PromptString','Select a figure to place cuts',...
    'SelectionMode','Single',...
    'ListSize',[300,160],...
    'ListString',AvaliableFigureStr);

% assign axes
if eq(Choise_hF_Plot,1)
    hF_Plot = figure;
    hAx_Plot = axes('Parent',hF_Plot);
else
    hF_Plot = hF_All(Choise_hF_Plot-1);
    hAx_Plot = findobj(hF_Plot,'Type','Axes');
end

% draw cuts
hold(hAx_Plot,'on')
    hL_Osc  =   plot(hAx_Plot,Osc_X,Osc_Y,'Color',Osc_C,'LineWidth',2);
    hIntL   = imline(hAx_imline,IntL_X,IntL_Y);
    assignin('base','hIntL',hIntL)
hold(hAx_Plot,'off')


hAx_Plot.XGrid = 'on';
hAx_Plot.YGrid = 'on';
setColor(hIntL,CursorColor)

hLs.hOsc    = hL_Osc;

addNewPositionCallback( hIntL,@(pos) SliceCallback(hAx_Plot,hLs,pos,h2D,Direction) );
        
switch Direction
    case 'X'        
        hAx_Plot.Tag = 'XCut';
        hAx_Plot.XLabel.String = XLabel;
        hAx_Plot.XLim = XLim;
        hAx_Plot.YLim = [-ZMax,ZMax];
        hAx_Plot.XMinorGrid = 'on';
    case 'Y'        
        hAx_Plot.Tag = 'YCut';
        hAx_Plot.YLabel.String = YLabel;
        hAx_Plot.XLim = YLim;
        hAx_Plot.YLim = [-ZMax,ZMax];
        hAx_Plot.YMinorGrid = 'on';
end

% Initialize Title
N_format   = '%6.2f';
CVlaue_Str = ['Current Point: (',num2str(X_Cut,N_format),', ',num2str(Y_Cut,N_format),')'];
hAx_Plot.Title.String = CVlaue_Str;

