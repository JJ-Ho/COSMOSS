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
defaultSaveFig     = 0;
defaultSavePath    = '~/Desktop/';
defaultFreqRange   = 1600:1800;
defaultNum_Contour = 20;
defaultPlotCursor  = 0;
defaultCMAP_Index  = 1;
defaultPlotNorm    = 0;

% add Optional inputs / Parameters
addOptional(INPUT,'SaveFig'    ,defaultSaveFig);
addOptional(INPUT,'SavePath'   ,defaultSavePath);
addOptional(INPUT,'FreqRange'  ,defaultFreqRange);
addOptional(INPUT,'Num_Contour',defaultNum_Contour);
addOptional(INPUT,'PlotCursor' ,defaultPlotCursor);
addOptional(INPUT,'CMAP_Index' ,defaultCMAP_Index);
addOptional(INPUT,'PlotNorm'   ,defaultPlotNorm);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
SaveFig     = INPUT.Results.SaveFig;
SavePath    = INPUT.Results.SavePath;
FreqRange   = INPUT.Results.FreqRange;
Num_Contour = INPUT.Results.Num_Contour;
PlotCursor  = INPUT.Results.PlotCursor;
CMAP_Index  = INPUT.Results.CMAP_Index;
PlotNorm    = INPUT.Results.PlotNorm;

%% Main
cla(hAx)
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
        set(hAx,'Ydir','normal')

    case 'None' % plot stick spectrum  w/ colormap
        spyXYZ(hAx,X,Y,Z_NC)    
        
        set(hAx,'Ydir','normal')
        hAx.Color = [0,0,0];
        
        colorbar
        StickC_Map = load('CoolWhite');
        colormap(StickC_Map.MAP) 
        Amp = max(abs(Z_NC(:)));
        caxis([-Amp,Amp])
        
    otherwise % plot convoluted spectrum
        contour(hAx,X,Y,Z,Num_Contour,'LineWidth',2)

        % Set colorbar
        colorbar
        CMAP = SelectColormap(CMAP_Index);
        colormap(CMAP)      
        Amp = max(abs(Z(:)));
        caxis([-Amp,Amp])
end

% Plot diagonal line
hold on; plot(hAx,X,Y); hold off

%% figure setting 
hF = hAx.Parent;
hF.Units = 'normalized'; % use normalized scale

hAx.DataAspectRatio = [1,1,1];
hAx.FontSize = 14;
hAx.XLabel.String = 'Probe (cm^{-1})';
hAx.YLabel.String = 'Pump (cm^{-1})';
hAx.XLim = [FreqRange(1),FreqRange(end)];
hAx.YLim = [FreqRange(1),FreqRange(end)];


FilesName       = CVL.FilesName;
FilesName_Reg   = regexprep(FilesName,'\_','\\_');
Coupling        = GUI_Inputs.CouplingType;
Coupling_Reg    = regexprep(Coupling,'\_','\\_');
Title_String{1} = [SpecType,' ',FilesName_Reg,', Coupling:',Coupling_Reg];
 
if PlotCursor
    Title_String{2} = '';
    % Call pointer
    S.fh = hF;
    S.ax = hAx;
    Pointer_N(S) % use normalized scale
end

title(Title_String,'FontSize',16);

%% Auto Save
if SaveFig
    timeStamp    = datetime('now','TimeZone','local');
    timeSamepStr = datestr(timeStamp,'yymmdd_HH-MM-SS');
    FigName      = [SpecType,'_',FilesName_Reg,'_',timeSamepStr];
    
    SaveFigures(hF,SavePath,FigName)
end