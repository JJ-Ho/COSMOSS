function  [SpectraGrid,Response] = TwoDIR_Main_Sparse(Structure,GUI_Inputs)
%% TwoDIR_Main(PDB_Data,GUI_Inputs)
%  
%   Given a initial stucture (pdb), this script will simulate its 2DIR
%   spectrum.
% 

% ------- Version log -----------------------------------------------------
% 
% Ver. 2.3  141014  Add "INCAng" as optional input.
% 
% Ver. 2.2  140723  Use input parser to remove GUI handles read-in. Now all
%                   of GUI read-in parts go into SFG_Acid.m layer.
% 
% Ver. 2.1  140717  Add Frequency read in function
% 
% Ver. 2.0  140716  Separate ploting into Convolution and making figure
% 
% Ver. 1.0  140716  Modified from TwoDSFG_AmideI_Main
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014

%% -Degub------------------------------------------
% PDB_Data = Structure;
% GUI_Inputs.debug = 1;
% %-Degub------------------------------------------

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

% Default values
defaultFreqRange    = 1650:1750;
defaultA_Pump1      = 90;
defaultA_Pump2      = 90;
defaultA_Probe      = 90;
defaultA_Sig2D      = 90;
defaultP_Pump1      = 0;
defaultP_Pump2      = 0;
defaultP_Probe      = 0;
defaultP_Sig2D      = 0;
defaultPCutOff      = 0;

addOptional(INPUT,'FreqRange'   ,defaultFreqRange);
addOptional(INPUT,'A_Pump1'     ,defaultA_Pump1);
addOptional(INPUT,'A_Pump2'     ,defaultA_Pump2);
addOptional(INPUT,'A_Probe'     ,defaultA_Probe);
addOptional(INPUT,'A_Sig2D'     ,defaultA_Sig2D);
addOptional(INPUT,'P_Pump1'     ,defaultP_Pump1);
addOptional(INPUT,'P_Pump2'     ,defaultP_Pump2);
addOptional(INPUT,'P_Probe'     ,defaultP_Probe);
addOptional(INPUT,'P_Sig2D'     ,defaultP_Sig2D);
addOptional(INPUT,'PCutOff'     ,defaultPCutOff);
         
parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
FreqRange    = INPUT.Results.FreqRange;
A_Pump1      = INPUT.Results.A_Pump1;
A_Pump2      = INPUT.Results.A_Pump2;
A_Probe      = INPUT.Results.A_Probe;
A_Sig2D      = INPUT.Results.A_Sig2D;
P_Pump1      = INPUT.Results.P_Pump1;
P_Pump2      = INPUT.Results.P_Pump2;
P_Probe      = INPUT.Results.P_Probe;
P_Sig2D      = INPUT.Results.P_Sig2D;
PCutOff      = INPUT.Results.PCutOff;

%% Call TwoExcitonH to calculate H,mu and alpha under exciton basis
H = ExcitonH(Structure,GUI_Inputs,'TwoEx');

Mu = MuAlphaGen(Structure,H,'Mode','Mu');

Ex_F1   = H.Sort_Ex_F1;
Ex_F2   = H.Sort_Ex_F2;
M_Ex_01 = Mu.M_Ex_01;
M_Ex_12 = Mu.M_Ex_12;

%% Decide what kinds of rod rotation average is and applied rotational 
% average on Response in molecular frame
R_Avg = LabFrameAvg('Isotropic',4); 

%% Jones Matrix convert XYZ to PS frame
% Turn degrees into radius (not work for BoxCard geometry yet)
A_Pump1 = A_Pump1/180*pi;
A_Pump2 = A_Pump2/180*pi;
A_Probe = A_Probe/180*pi;
A_Sig2D = A_Sig2D/180*pi;

J = JonesTrans4(A_Sig2D,A_Probe,A_Pump2,A_Pump1);

%% E part, Plarization of each incident beams
% Polarization Angles of incident beams
P_Pump1  = P_Pump1/180*pi;
P_Pump2  = P_Pump2/180*pi;
P_Probe  = P_Probe/180*pi;
P_Sig2D  = P_Sig2D/180*pi;

E = EPolar4(P_Sig2D,P_Probe,P_Pump2,P_Pump1);

%% Generate Feynman pathway for 2DSFG
EJR = E*J*R_Avg;

[Grid,Freq,Int,Index,CutOff] = Feynman_2DIR_Vec_Sparse(PCutOff,...
                                                       FreqRange,...
                                                       EJR,...
                                                       Ex_F1,...
                                                       Ex_F2,...
                                                       M_Ex_01,...
                                                       M_Ex_12);

%% Group up outputs
Response.H = H;
Response.Mu = Mu;

Response.Freq  = Freq;
Response.Int   = Int;
Response.Index = Index;
Response.CutOff = CutOff;

Response.EJR = EJR;
Response.SpecType = '2DIR';

%% Binning of spectra
SpectraGrid.SpecAccuR1  = Grid.R1;
SpectraGrid.SpecAccuR2  = Grid.R2;
SpectraGrid.SpecAccuR3  = Grid.R3;
SpectraGrid.SpecAccuNR1 = Grid.NR1;
SpectraGrid.SpecAccuNR2 = Grid.NR2;
SpectraGrid.SpecAccuNR3 = Grid.NR3;

SpectraGrid.Rephasing    = Grid.R1  + Grid.R2  - Grid.R3;
SpectraGrid.NonRephasing = Grid.NR1 + Grid.NR2 - Grid.NR3;

