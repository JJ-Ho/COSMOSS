function hF = Plot2D(hAx,CVL,GUI_Inputs,SpecType)
% 
% This function plot 2DIR and other information

% ------- Version log -----------------------------------------------------
%  
% Ver. 1.2  161013  Add handle of figure for dynamic updates
% 
% Ver. 1.1  141014  Use Pointer_N instead of Pointer for normailized unit
% 
% Ver. 1.0  140723  Isolated from "TwoDIR_Main.m"
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014-2016

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultSaveFig        = 0;
defaultSavePath       = '~/Desktop/';
defaultexistFig       = 0;
defaulthFig           = '';
defaultPlotNorm_2D    = 0;
defaultPlotCursor_2D  = 0;
defaultNum_Contour_2D = 20;
defaultCMAP_Index_2D  = 1;
defaultPathway_2D     = 'All';

% add Optional inputs / Parameters
addOptional(INPUT,'SaveFig'       ,defaultSaveFig);
addOptional(INPUT,'SavePath'      ,defaultSavePath);
addOptional(INPUT,'existFig'      ,defaultexistFig);
addOptional(INPUT,'hFig'          ,defaulthFig);
addOptional(INPUT,'PlotNorm_2D'   ,defaultPlotNorm_2D);
addOptional(INPUT,'PlotCursor_2D' ,defaultPlotCursor_2D);
addOptional(INPUT,'Num_Contour_2D',defaultNum_Contour_2D);
addOptional(INPUT,'CMAP_Index_2D' ,defaultCMAP_Index_2D);
addOptional(INPUT,'Pathway_2D'  ,defaultPathway_2D  );

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
SaveFig     = INPUT.Results.SaveFig;
SavePath    = INPUT.Results.SavePath;
existFig    = INPUT.Results.existFig;
hFig        = INPUT.Results.hFig;
PlotCursor  = INPUT.Results.PlotCursor_2D;
PlotNorm    = INPUT.Results.PlotNorm_2D;
Num_Contour = INPUT.Results.Num_Contour_2D;
CMAP_Index  = INPUT.Results.CMAP_Index_2D;
Pathway     = INPUT.Results.Pathway_2D;

%% Main
FreqRange = CVL.FreqRange;

if existFig
    hAx = findobj(evalin('base',hFig),'Type','Axes');
end

hold(hAx,'on')
X = FreqRange;
Y = FreqRange;

NC = CVL.selected_No_Conv;
C  = real(CVL.selected);

if PlotNorm
    Z_NC = NC ./ max(abs(NC(:)));
    Z    =  C ./ max(abs( C(:)));
else
    Z_NC = NC;
    Z    =  C;
end

switch CVL.Lineshape
    case 'Spy' % plot stick spectrum
        spyXY(hAx,X,Y,Z_NC)  
        
        C_MAP = SelectColormap('Sign');
        C_Amp = 1;
        DiagColor = 'k';
        hAx.YDir  = 'normal';

    case 'No Conv' % plot stick spectrum  w/ colormap
        spyXYZ(hAx,X,Y,Z_NC)    

        C_MAP = SelectColormap('CoolWhite');hAx.Color = [0,0,0];
        %C_MAP = SelectColormap('JetBlackWide') ;hAx.Color = [1,1,1];
        C_Amp = max(abs(Z_NC(:)));
        DiagColor = 'g';
        hAx.YDir  = 'normal';

    otherwise % plot convoluted spectrum
        [~,hC] = contourf(hAx,X,Y,Z,Num_Contour,'LineWidth',2);
        
        hC.LineStyle = 'none';
        C_MAP = SelectColormap(CMAP_Index);      
        C_Amp = max(abs(Z(:)));
        DiagColor = 'k';
        
        Fig_ContourLine = 1;
        if Fig_ContourLine
            % copy the contour and make it into line only
            hCL = copyobj(hC,hAx);
            hCL.Fill = 'off';
            hCL.LineStyle = '-';
            hCL.LineWidth = 0.5;
            hCL.LineColor = [0,0,0]; 
            NC_half = floor(Num_Contour/2);
            CLevelList = linspace(-1,1,NC_half)*C_Amp;
            %Fig_CL_D_Ind = [floor(NC_half/2):floor(NC_half/2)+1];
            Fig_CL_D_Ind = floor(NC_half/2)+1;
            CLevelList(Fig_CL_D_Ind)=[]; % Deal with deleting zero contour line
            hCL.LevelList = CLevelList;
        end
end
colorbar(hAx)
colormap(hAx,C_MAP)
caxis(hAx,[-C_Amp,C_Amp])

% Plot diagonal line
plot(hAx,X,Y,'Color',DiagColor,'LineStyle','--')
hold(hAx,'off')

%% figure setting 
hAx.DataAspectRatio = [1,1,1];
hAx.FontSize = 14;
hAx.XLabel.String = 'Probe (cm^{-1})';
hAx.YLabel.String = 'Pump (cm^{-1})';
hAx.XLim = [FreqRange(1),FreqRange(end)];
hAx.YLim = [FreqRange(1),FreqRange(end)];

FilesName     = CVL.FilesName;
FilesName_Reg = regexprep(FilesName,'\_','\\_');
% Coupling      = GUI_Inputs.CouplingType;
% Coupling_Reg  = regexprep(Coupling,'\_','\\_');
%Title_String = [SpecType,' ',FilesName_Reg,', Coupling:',Coupling_Reg];

Title_String = [SpecType,' ',FilesName_Reg,', Pathway: ',Pathway];
hAx.Title.String   = Title_String;
hAx.Title.FontSize = 16;

if PlotCursor
    hF = hAx.Parent;
    hF.Units = 'normalized'; % use normalized scale
    S.hF = hF;
    S.hAx = hAx;
    Pointer_T(S);
end

%% Auto Save
if SaveFig
    timeStamp    = datetime('now','TimeZone','local');
    timeSamepStr = datestr(timeStamp,'yymmdd_HH-MM-SS');
    FigName      = [SpecType,'_',timeSamepStr];
    
    hF = hAx.Parent;
    SaveFigures(hF,SavePath,FigName)
end

