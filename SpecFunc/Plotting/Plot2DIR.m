function hF = Plot2DIR(CVL,GUI_Inputs)
% 
% This function plot 2DIR and other information

% ------- Version log -----------------------------------------------------
%  
% Ver. 1.1  141014  Use Pointer_N instead of Pointer for normailized unit
% 
% Ver. 1.0  140723  Isolated from "TwoDIR_Main.m"
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014

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

% add Optional inputs / Parameters
addOptional(INPUT,'FreqRange'  ,defaultFreqRange);
addOptional(INPUT,'Num_Contour',defaultNum_Contour);
addOptional(INPUT,'PlotCursor' ,defaultPlotCursor);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
FreqRange   = INPUT.Results.FreqRange;
Num_Contour = INPUT.Results.Num_Contour;
PlotCursor  = INPUT.Results.PlotCursor;

%% Main
hF = figure;

CVLRS = real(CVL.selected);
contour(FreqRange,FreqRange,CVLRS,Num_Contour,'LineWidth',2)

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

% Set colorbar
colorbar
colormap('jet')
Amp = max(abs(caxis));
caxis([-Amp Amp])

if PlotCursor
    % Call pointer
    S.fh = hF;
    S.ax = hAx;
    Pointer_N(S) % use normalized scale
else
    FilesName     = CVL.FilesName;
    FilesName_Reg = regexprep(FilesName,'\_','\\_');
    Coupling      = GUI_Inputs.Coupling;
    Coupling_Reg  = regexprep(Coupling,'\_','\\_');
    Title_String  = ['2DIR ',FilesName_Reg,', Coupling:',Coupling_Reg];
    title(Title_String,'FontSize',16); 
end
