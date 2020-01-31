function Output = FTIR_Main(Structure,GUI_Inputs)
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
defaultFreqRange    = 1650:1750;

addOptional(INPUT,'FreqRange'   ,defaultFreqRange);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
FreqRange    = INPUT.Results.FreqRange;

%% Main
H = ExcitonH(Structure,GUI_Inputs,'OneEx');

Mu = MuAlphaGen(Structure,H,'Mode','Mu');

Ex_F1   = Mu.Sort_Ex_F1;
M_Ex_01 = Mu.M_Ex_01;

Response = sum(M_Ex_01.^2,2); % E-field of FTIR signal is mu^2 base on feynmann diagram! 

%% Bin signal
AccuGrid = Bin1D(Ex_F1,Response,FreqRange);

%% Output
Output.FilesName    = Structure.FilesName;
Output.SpecType     = 'FTIR';
Output.Response1D   = AccuGrid;
Output.freq_OneD    = FreqRange;
Output.H            = H;
Output.Mu           = Mu;




