function  [SpectraGrid,Response] = TwoDIR_Main(PDB_Data,GUI_Inputs)
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
defaultCouplingType = 'TDC'; 
% expectedCoupling   = {'NN_Mix_TDC','TDC','Cho_PB','Cho_APB'};
defaultBeta_NN      = 0.8; % 0.8 cm-1 according to Lauren's PNAS paper (doi/10.1073/pnas.1117704109); that originate from Min Cho's paper (doi:10.1063/1.1997151)
defaultA_Pump1      = 90;
defaultA_Pump2      = 90;
defaultA_Probe      = 90;
defaultA_Sig2D      = 90;
defaultP_Pump1      = 0;
defaultP_Pump2      = 0;
defaultP_Probe      = 0;
defaultP_Sig2D      = 0;


addOptional(INPUT,'FreqRange'   ,defaultFreqRange);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);
addOptional(INPUT,'A_Pump1'     ,defaultA_Pump1);
addOptional(INPUT,'A_Pump2'     ,defaultA_Pump2);
addOptional(INPUT,'A_Probe'     ,defaultA_Probe);
addOptional(INPUT,'A_Sig2D'     ,defaultA_Sig2D);
addOptional(INPUT,'P_Pump1'     ,defaultP_Pump1);
addOptional(INPUT,'P_Pump2'     ,defaultP_Pump2);
addOptional(INPUT,'P_Probe'     ,defaultP_Probe);
addOptional(INPUT,'P_Sig2D'     ,defaultP_Sig2D);
addOptional(INPUT,'CouplingType',defaultCouplingType);
% addParamValue(INPUT,'Coupling',defaultCoupling,...
%                  @(x) any(validatestring(x,expectedCoupling)));
         
parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
FreqRange    = INPUT.Results.FreqRange;
CouplingType = INPUT.Results.CouplingType;
Beta_NN      = INPUT.Results.Beta_NN;
A_Pump1      = INPUT.Results.A_Pump1;
A_Pump2      = INPUT.Results.A_Pump2;
A_Probe      = INPUT.Results.A_Probe;
A_Sig2D      = INPUT.Results.A_Sig2D;
P_Pump1      = INPUT.Results.P_Pump1;
P_Pump2      = INPUT.Results.P_Pump2;
P_Probe      = INPUT.Results.P_Probe;
P_Sig2D      = INPUT.Results.P_Sig2D;

%% Call TwoExcitonH to calculate H,mu and alpha under exciton basis
H = ExcitonH(PDB_Data,...
             'ExMode'  ,'TwoEx',...
             'CouplingType',CouplingType,...
             'Beta_NN' ,Beta_NN);

% Mu = MuAlphaGen_full_M(PDB_Data,H,'Mode','Mu');
% Sort_Ex_Freq = H.Sort_Ex_Freq;
% Mu_Ex        = Mu.Trans_Ex;

Mu = MuAlphaGen(PDB_Data,H,'Mode','Mu');

Ex_F1   = H.Sort_Ex_F1;
Ex_F2   = H.Sort_Ex_F2;
M_Ex_01 = Mu.M_Ex_01;
M_Ex_12 = Mu.M_Ex_12;

%% Generate Feynman pathway for 2DSFG
Num_Modes = PDB_Data.Num_Modes;
% Response = Feynman_2DIR_kron(Num_Modes,Sort_Ex_Freq,Mu_Ex);
% Response = Feynman_2DIR_Vec_Full_M(Num_Modes,Sort_Ex_Freq,Mu_Ex);
[Freq,Beta,Index] = Feynman_2DIR_Vec(Num_Modes,Ex_F1,Ex_F2,M_Ex_01,M_Ex_12);

%% Decide what kinds of rod rotation average is and applied rotational 
% average on Response in molecular frame
R_Avg = LabFrameAvg('Isotropic',4); 

AR1  = R_Avg*Beta.R1 ;
AR2  = R_Avg*Beta.R2 ;
AR3  = R_Avg*Beta.R3 ;
NAR1 = R_Avg*Beta.NR1;
NAR2 = R_Avg*Beta.NR2;
NAR3 = R_Avg*Beta.NR3;

%% Jones Matrix convert XYZ to PS frame
% Turn degrees into radius (not work for BoxCard geometry yet)
A_Pump1 = A_Pump1/180*pi;
A_Pump2 = A_Pump2/180*pi;
A_Probe = A_Probe/180*pi;
A_Sig2D = A_Sig2D/180*pi;

J = JonesTrans4(A_Sig2D,A_Probe,A_Pump2,A_Pump1);

JAR1  = J*AR1;
JAR2  = J*AR2;
JAR3  = J*AR3;
JNAR1 = J*NAR1;
JNAR2 = J*NAR2;
JNAR3 = J*NAR3;

%% E part, Plarization of each incident beams
% Polarization Angles of incident beams
P_Pump1  = P_Pump1/180*pi;
P_Pump2  = P_Pump2/180*pi;
P_Probe  = P_Probe/180*pi;
P_Sig2D  = P_Sig2D/180*pi;

E = EPolar4(P_Sig2D,P_Probe,P_Pump2,P_Pump1);

EJAR1  = E*JAR1;
EJAR2  = E*JAR2;
EJAR3  = E*JAR3;
EJNAR1 = E*JNAR1;
EJNAR2 = E*JNAR2;
EJNAR3 = E*JNAR3;

%% Output
Response.H  = H;
Response.Mu = Mu;

Response.Freq  = Freq;
Response.Beta  = Beta;
Response.Index = Index;

Response.RBeta.R_Avg = R_Avg;
Response.RBeta.R1    = AR1;
Response.RBeta.R2    = AR2;
Response.RBeta.R3    = AR3;
Response.RBeta.NR1   = NAR1;
Response.RBeta.NR2   = NAR2;
Response.RBeta.NR3   = NAR3;

Response.JRBeta.J    = J;
Response.JRBeta.R1   = JAR1;
Response.JRBeta.R2   = JAR2;
Response.JRBeta.R3   = JAR3;
Response.JRBeta.NR1  = JNAR1;
Response.JRBeta.NR2  = JNAR2;
Response.JRBeta.NR3  = JNAR3;

Response.EJRBeta.E   = E;
Response.EJRBeta.R1  = EJAR1;
Response.EJRBeta.R2  = EJAR2;
Response.EJRBeta.R3  = EJAR3;
Response.EJRBeta.NR1 = EJNAR1;
Response.EJRBeta.NR2 = EJNAR2;
Response.EJRBeta.NR3 = EJNAR3;

%% Binning of spectra
S = Response.EJRBeta;
SpectraGrid  = Bin2D(S,Freq,FreqRange);
