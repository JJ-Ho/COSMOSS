function hF = Plot2DIR(CVL,H,Mu_Ex,GUI_Inputs)
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
defaultPlotCursor  = 0;

% add Optional inputs / Parameters
addOptional(INPUT,'FreqRange'  ,defaultFreqRange);
addOptional(INPUT,'PlotCursor' ,defaultPlotCursor);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
FreqRange   = INPUT.Results.FreqRange;
PlotCursor  = INPUT.Results.PlotCursor;

%% Main
hF = figure('Position',[600,200,500,450],...
           'units','pixel');

% Generate Mu value table and eigen mode frequencies
Mod_Ine_All = 2:3; % for two coupled oscillators

F = H.Sort_Ex_Freq(Mod_Ine_All)';


Mu_Vec = squeeze(Mu_Ex(1,Mod_Ine_All,:));
Mu_Vec_L = sqrt(sum(Mu_Vec.^2,2))';
T = [ Mu_Vec_L.^2;
      Mu_Vec_L.^4;
      Mu_Vec(:,1)';
      Mu_Vec(:,2)';
      Mu_Vec(:,3)'];

Beta = ones(1,2).*H.Beta(1,2);
  
dat = [F;Beta;T]'; 
rnames = {'S1','S2'};
cnames = {'Freq','Beta','IR_Int','2DIR_Int','Mu_X','Mu_Y','Mu_Z'};

t = uitable('Parent',hF,...
            'Data',dat,...
            'ColumnFormat',{'bank','bank','bank','bank','bank','bank','bank'},...
            'ColumnWidth' ,{60,50,60,80,40,40,40},...
            'ColumnName',cnames,... 
            'RowName',rnames,...
            'units','pix',...
            'fontsize',13,...
            'Position',[40,10,420,60]);
        
% Make 2DSFG plot
hAx = axes('units','pix','Position',[100,80,350,350]);
% h1 = axes('units','pix','Position',[100,100,300,300]);
set(gcf,'CurrentAxes',hAx)

% Num_Contour = 20;
% [C,Ch] = contour(FreqRange,FreqRange,real(CVL.selected),Num_Contour,'LineWidth',2);
% shading flat
% set(gca,'DataAspectRatio',[1 1 1])

Z = -1:(2/40):1;
% V = Z;
V = [Z(1:20) Z(22:41)];

% contour(FreqRange,FreqRange,-1.*real(CVL.selected)./max(max(abs(real(CVL.selected)))),V,'LineWidth',1.5)
contour(FreqRange,FreqRange,-1.*real(CVL.selected),40,'LineWidth',1.5)
% caxis([-1 1]);
shading flat
set(gca,'DataAspectRatio',[1 1 1])

% Plot diagonal line
hold on; plot(FreqRange,FreqRange); hold off

%% figure setting 
hAx.FontSize = 14;
hAx.XLabel.String = 'Probe (cm^{-1})';
hAx.YLabel.String = 'Pump (cm^{-1})';

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
    Title = ['2DIR'];
    title(Title,'FontSize',16);    
end

% figure; plot(FreqRange,diag(real(CVL.selected)))