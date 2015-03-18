function Output = PlotFTIR(PDB_Data,varargin)
%% PlotFTIR
%  
% Given one exciton alpha, mu matrix and respective one exciton state
% energy (cm-1), this function can generate 1DSFG plots with different
% polarization combinations {xx,yy,zz,xy,yz,xz}*{x,y,z} = 18 plots.
% 
% 
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.4  140922  Add Output part;
%                   Add Inputparser
% 
% Ver. 1.3  140717  Add Frequency axis GUI read in part
% 
% Ver. 1.2  140616  Add CouplingOpt
% 
% Ver. 1.1  140605  Isolate from PlotOneDSFG.m
% 
% Ver. 1.0  130729  Isolated from TwoDSFG_Simulation.
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

%% Debug
% PDB_Data = GetAmideI;
% handles.PDB_Data = PDB_Data;
%% Inputs parser
INPUT = inputParser;

% Default values
defaultLabel_Index = 'Non';
defaultLabel_Freq  = 1700;
defaultPlotStick   = 1;
defaultCoupling    = 1;
defaultBeta_NN     = 0.8;
defaultF_Min       = 1600;
defaultF_Max       = 1800;
defaultLineWidth   = 5;

% add Optional inputs / Parameters
addOptional(INPUT,'Label_Index',defaultLabel_Index);
addOptional(INPUT,'Label_Freq' ,defaultLabel_Freq);
addOptional(INPUT,'PlotStick'  ,defaultPlotStick);
addOptional(INPUT,'Coupling'   ,defaultCoupling);
addOptional(INPUT,'Beta_NN'    ,defaultBeta_NN);
addOptional(INPUT,'F_Min'      ,defaultF_Min);
addOptional(INPUT,'F_Max'      ,defaultF_Max);
addOptional(INPUT,'LineWidth'  ,defaultLineWidth);

parse(INPUT,varargin{:});

% Re-assign variable names
Label_Index = INPUT.Results.Label_Index;
Label_Freq  = INPUT.Results.Label_Freq;
PlotStick   = INPUT.Results.PlotStick;
Coupling    = INPUT.Results.Coupling;
Beta_NN     = INPUT.Results.Beta_NN;
F_Min       = INPUT.Results.F_Min;
F_Max       = INPUT.Results.F_Max;
LineWidth   = INPUT.Results.LineWidth;

%% Read GUI inputs
% PDB_Data = handles.PDB_Data;
% Label_Index = str2num(get(handles.E_LIndex,'String'));
% Label_Freq  = str2double(get(handles.E_LFreq,'String'));

% PlotStick = get(handles.PlotStick,'Value');
% Coupling = get(handles.Coupling,'Value');

% F_Min = str2double(get(handles.X_Min,'String'));
% F_Max = str2double(get(handles.X_Max,'String'));

% LineWidth=str2double(get(handles.LineWidth,'String'));
%% Main

Num_Modes = PDB_Data.Num_Modes;

% isotope labled frequencies
if ~ischar(Label_Index)
    PDB_Data.freq(Label_Index) = Label_Freq.*ones(size(Label_Index));
end

switch Coupling
    case 1
        Coupling = 'TDC';
    case 2
        Coupling = 'Zero';
    case 3
        Coupling = 'NN';
    case 4 
        Coupling = 'NN_Mix_TDC';
end

H = ExcitonH(PDB_Data,'ExMode','OneEx','Coupling',Coupling,'Beta_NN',Beta_NN);

Freq01Ex = H.Sort_Ex_Freq;

mu    = MuAlphaGen(PDB_Data,H,'Mode','Mu');
mu01Ex = mu.Trans_Ex;

%Intensity for linear spec
% IntM = sqrt(sum(mu01Ex.^2,3));
IntM = sum(mu01Ex.^2,3);

% Test plot of 1D FTIR
mu_OneD = IntM(2:Num_Modes+1,1);
freq_OneD = Freq01Ex(2:Num_Modes+1);

f = figure; hold on
set(f,'Unit','normalized') % use normalized scale

if eq(PlotStick,1)
%     plot(freq_OneD,mu_OneD,'rx')
    line([freq_OneD';freq_OneD'],[zeros(1,Num_Modes);mu_OneD'])
end
    
% Get Frequency axis range
spec_range = F_Min:F_Max;

spec_array1 = bsxfun(@times,ones(Num_Modes,length(spec_range)),spec_range);
spec_array2 = bsxfun(@minus,spec_array1,freq_OneD);

Gaussian = bsxfun(@times,exp(-(spec_array2.^2)./(LineWidth^2)),mu_OneD);
Gaussian_Toatl = sum(Gaussian,1);

plot(spec_range,Gaussian_Toatl,'-')
set(gca,'XLim',[spec_range(1),spec_range(end)])

hold off

% % Call pointer
% S.fh = f;
% S.ax = get(f,'CurrentAxes');
% Pointer_N(S) % use normalized scale

% integrate the curve area
Area = trapz(spec_range,Gaussian_Toatl);
StickSum = sum(mu_OneD);

Title = [num2str(Area), ',  ', num2str(StickSum)];
title(Title,'FontSize',16);


%% Output
Output.H  = H;
Output.Mu = mu;
