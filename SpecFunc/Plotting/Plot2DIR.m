function hF = Plot2DIR(hF,CVL,GUI_Inputs)
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
defaultFreqRange   = 1600:1800;
defaultNum_Contour = 20;
defaultPlotCursor  = 0;
defaultCMAP_Index  = 1;

% add Optional inputs / Parameters
addOptional(INPUT,'FreqRange'  ,defaultFreqRange);
addOptional(INPUT,'Num_Contour',defaultNum_Contour);
addOptional(INPUT,'PlotCursor' ,defaultPlotCursor);
addOptional(INPUT,'CMAP_Index' ,defaultCMAP_Index);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
FreqRange   = INPUT.Results.FreqRange;
Num_Contour = INPUT.Results.Num_Contour;
PlotCursor  = INPUT.Results.PlotCursor;
CMAP_Index  = INPUT.Results.CMAP_Index;

%% Main
hAx = findobj(hF,'type','axes');
cla(hAx)

if strcmp(CVL.Lineshape,'None')
    % plot stick spectrum
    imagesc(FreqRange,FreqRange,CVL.selected_No_Conv)
    set(gca,'Ydir','normal')
    
    % Set colorbar
    colorbar
    StickC_Map = load('CoolBlack');
    colormap(StickC_Map.MAP) 
    Amp = max(abs(CVL.selected_No_Conv(:)));
    caxis([-Amp,Amp])
else
    % plot convoluted spectrum
    CVLRS = real(CVL.selected);
    contour(FreqRange,FreqRange,CVLRS,Num_Contour,'LineWidth',2)
    % Normalization
    % CVLRSN = CVLRS ./max(abs(CVLRS(:)));
    % contour(GUI_Inputs.FreqRange,GUI_Inputs.FreqRange,CVLRSN,GUI_Inputs.Num_Contour,'LineWidth',2)
    
    % Set colorbar
    colorbar
    CMAP = SelectColormap(CMAP_Index);
    colormap(CMAP)      
    Amp = max(abs(CVLRS(:)));
    caxis([-Amp,Amp])
end

% Plot diagonal line
hold on; plot(FreqRange,FreqRange); hold off

%% figure setting 
hF.Units = 'normalized'; % use normalized scale
hAx = hF.CurrentAxes;
hAx.DataAspectRatio = [1,1,1];
hAx.FontSize = 14;
hAx.XLabel.String = 'Probe (cm^{-1})';
hAx.YLabel.String = 'Pump (cm^{-1})';

shading flat

if PlotCursor
    % Call pointer
    S.fh = hF;
    S.ax = hAx;
    Pointer_N(S) % use normalized scale
else
    FilesName     = CVL.FilesName;
    FilesName_Reg = regexprep(FilesName,'\_','\\_');
    Coupling      = GUI_Inputs.CouplingType;
    Coupling_Reg  = regexprep(Coupling,'\_','\\_');
    Title_String  = ['2DIR ',FilesName_Reg,', Coupling:',Coupling_Reg];
    title(Title_String,'FontSize',16); 
end
