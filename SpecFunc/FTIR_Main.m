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
defaultFreqRange_1D    = 1650:1750;
defaultAvg_Rot      = 1;
defaultAvg_Mirror   = 1;

addOptional(INPUT,'FreqRange_1D',defaultFreqRange_1D);
addOptional(INPUT,'Avg_Rot'     ,defaultAvg_Rot);
addOptional(INPUT,'Avg_Mirror'  ,defaultAvg_Mirror);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
FreqRange    = INPUT.Results.FreqRange_1D;
Avg_Rot      = INPUT.Results.Avg_Rot;
Avg_Mirror   = INPUT.Results.Avg_Mirror;

%% Main
H  = H_handler(Structure,GUI_Inputs,'OneEx');
Mu = TrMoment(Structure,'OneEx','Mu',H);

Ex_F1   = H.Sort_Ex_F1;
M_Ex_01 = Mu.M_Ex_01'; %=> [3*N]

%% Generate Molecular frame SFG Responses
% vectorized version
[M_Ind1,M_Ind2] = ndgrid(1:3,1:3);
ResLF = M_Ex_01(M_Ind1(:),:).* M_Ex_01(M_Ind2(:),:); %=> [9,N]
ResLF_Int = sqrt(sum(ResLF.^2,1)); %=> [1,N]

%% Decide what kinds of ensemble average
N_Interactions = 2; % for FTIR
[R_Avg,Mirror_Mask,~,~] = LabFrameAvg(Avg_Rot,Avg_Mirror,N_Interactions);

ResLF_Avg = bsxfun(@times,R_Avg*ResLF,Mirror_Mask); %=> [9,N]

%% Jones Matrix convert XYZ to PS frame
% 
% JLabFrame = [freq, pp, ps, sp, ss]
% 
% Turn degrees into radius
% Normal incident and collection:
A_IR_IN  = 0;
A_IR_OUT = 0;

A_IR_IN  =  A_IR_IN/180*pi;
A_IR_OUT = A_IR_OUT/180*pi;

J = JonesTrans2(A_IR_OUT,A_IR_IN); % => [9,4]

J_ResLF_Avg = J * ResLF_Avg; % => [4,N]

%% E part, Plarization of each incident beams
% Polarization Angles of incident beams, 0 = P, 90 = S
P_IR = 0;
P_IR = P_IR/180*pi;

E = EPolar2(P_IR,P_IR); % => [1,4]

E_J_ResLF_Avg = E * J_ResLF_Avg; % => [1,N]

%% Bin signal
AccuGrid    = Bin1D(Ex_F1,E_J_ResLF_Avg,FreqRange);
AccuGridInt = Bin1D(Ex_F1,    ResLF_Int,FreqRange);

%% Output
Output.FilesName    = Structure.FilesName;
Output.SpecType     = 'FTIR';
Output.Response1D   = AccuGrid;
Output.Res_Int      = AccuGridInt;
Output.freq_OneD    = FreqRange;
Output.Mu           = Mu;
Output.H            = H;
Output.MolFrame     = ResLF;
Output.R_Avg        = R_Avg;
Output.LabFrame     = ResLF_Avg;
Output.Jones        = J;
Output.JLabFrame    = J_ResLF_Avg;
Output.E            = E;
Output.EJLabFrame   = E_J_ResLF_Avg;

