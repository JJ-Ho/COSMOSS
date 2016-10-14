function Output = PlotFTIR(PDB_Data,GUI_Inputs)
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
% PDB_Data = GetAcid;
% handles.PDB_Data = PDB_Data;
% 
% GUI_Inputs.debug = 1;
% GUI_Inputs.Coupling = 'NN';

%% Inputs parser
% Turn Output from Read GUI to cell array
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

% Default values
defaultCouplingType = 'TDC';
defaultBeta_NN      = 0.8;
defaultFreqRange    = 1650:1750;


% add Optional inputs / Parameters
addOptional(INPUT,'CouplingType',defaultCouplingType);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);
addOptional(INPUT,'FreqRange'   ,defaultFreqRange);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
CouplingType = INPUT.Results.CouplingType;
Beta_NN      = INPUT.Results.Beta_NN;
FreqRange    = INPUT.Results.FreqRange;

%% Main
Num_Modes = PDB_Data.Num_Modes;

H = ExcitonH(PDB_Data,'ExMode','OneEx','CouplingType',CouplingType,'Beta_NN',Beta_NN);

mu = MuAlphaGen(PDB_Data,H,'Mode','Mu');

muEx = mu.Trans_Ex;
muEx_Vec = muEx(2:Num_Modes+1,1,:);

Freq01Ex = H.Sort_Ex_Freq;
freq_OneD = Freq01Ex(2:Num_Modes+1);
freq_OneD = round(freq_OneD); % bin the fre to 1cm^-1

mu2_OneD = sum(muEx_Vec.^2,3); % E-field of FTIR signal is mu^2 base on feynmann duagram! 

%% Bin signal
AccuGrid = Bin1D(freq_OneD,mu2_OneD,FreqRange);


%% Output
Output.Num_Modes    = Num_Modes;
Output.Response1D   = AccuGrid;
Output.freq_OneD    = FreqRange;
Output.SpecType     = 'FTIR';
Output.FilesName    = PDB_Data.FilesName;
Output.CouplingType = CouplingType;
Output.H            = H;
Output.Mu           = mu;
