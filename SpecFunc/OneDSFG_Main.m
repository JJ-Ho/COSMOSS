function  OneDSFG = OneDSFG_Main(Structure,GUI_Inputs)
%% TwoDSFG_AmideI_Main

% ------- Version log -----------------------------------------------------
% 
% Ver. 2.0  161108  Vectorized version, separate frequency from calculated
%                   response matrix
% 
% Ver. 1.3  140922  Add input parser back
% 
% Ver. 1.2  140607  Delete input parser
% 
% Ver. 1.1  140603  Integrate isotope labeling from GetAmideI.m
% 
% Ver. 1.0  140603  Modified from TwoDSFG_AmideI for OneDSFG_AmideI GUI
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014

%% -Degub------------------------------------------
% clear all
% Structure_Data = GetAcid;
% GUI_Inputs.debug = 1;
% -Degub------------------------------------------

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

% Default values
defaultCouplingType = 'TDC';
defaultAvg_Rot      = 1;
defaultAvg_Mirror   = 1;
defaultA_IR         = 90;
defaultA_Vis1D      = 90;
defaultA_Sig1D      = 90;
defaultP_IR         = 0;
defaultP_Vis1D      = 0;
defaultP_Sig1D      = 0;
% defaultBeta_NN      = 0.8;
defaultFreqRange    = 1650:1750;

% add Optional inputs / Parameters
addOptional(INPUT,'CouplingType',defaultCouplingType);
addOptional(INPUT,'Avg_Rot'     ,defaultAvg_Rot);
addOptional(INPUT,'Avg_Mirror'  ,defaultAvg_Mirror);
addOptional(INPUT,'A_IR'        ,defaultA_IR);
addOptional(INPUT,'A_Vis1D'     ,defaultA_Vis1D);
addOptional(INPUT,'A_Sig1D'     ,defaultA_Sig1D);
addOptional(INPUT,'P_IR'        ,defaultP_IR);
addOptional(INPUT,'P_Vis1D'     ,defaultP_Vis1D);
addOptional(INPUT,'P_Sig1D'     ,defaultP_Sig1D);
% addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);
addOptional(INPUT,'FreqRange'   ,defaultFreqRange);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names

% CouplingType = INPUT.Results.CouplingType;
Avg_Rot      = INPUT.Results.Avg_Rot;
Avg_Mirror   = INPUT.Results.Avg_Mirror;
A_IR         = INPUT.Results.A_IR;
A_Vis1D      = INPUT.Results.A_Vis1D;
A_Sig1D      = INPUT.Results.A_Sig1D;
P_IR         = INPUT.Results.P_IR;
P_Vis1D      = INPUT.Results.P_Vis1D;
P_Sig1D      = INPUT.Results.P_Sig1D;
% Beta_NN      = INPUT.Results.Beta_NN;
FreqRange    = INPUT.Results.FreqRange;

%% Call OneExcitonH to calculate H,mu and alpha under exciton basis
% H = ExcitonH(SData,'ExMode','OneEx','CouplingType',CouplingType,'Beta_NN',Beta_NN);
H = ExcitonH(Structure,GUI_Inputs,'OneEx');

Mu    = MuAlphaGen(Structure,H,'Mode','Mu');
Alpha = MuAlphaGen(Structure,H,'Mode','Alpha');

Ex_F1   = H.Sort_Ex_F1;
M_Ex_01 = Mu.M_Ex_01';    %=> [3*N]
A_Ex_01 = Alpha.M_Ex_01'; %=> [9*N]

%% Generate Molecular frame SFG Responses
% vectorized version
[M_Ind,A_Ind] = ndgrid(1:3,1:9);
ResLF = A_Ex_01(A_Ind(:),:).* M_Ex_01(M_Ind(:),:); %=> [27*N]

%% Decide what kinds of ensemble average
N_Interactions = 3; % for SFG
[R_Avg,Mirror_Mask,~,~] = LabFrameAvg(Avg_Rot,Avg_Mirror,N_Interactions);

ResLF_Avg = bsxfun(@times,R_Avg*ResLF,Mirror_Mask); %=> [27*N]

%% Jones Matrix convert XYZ to PS frame
% 
% JLabFrame = [freq, ppp, pps, psp, pss, spp, sps, ssp, sss]
% 
% Turn degrees into radius
A_IR    =    A_IR/180*pi;
A_Vis1D = A_Vis1D/180*pi;
A_Sig1D = A_Sig1D/180*pi;

J = JonesRef3(A_Sig1D,A_Vis1D,A_IR); % => [27,8]

J_ResLF_Avg = J * ResLF_Avg; % => [8,N]

%% E part, Plarization of each incident beams
% Polarization Angles of incident beams, 0 = P, 90 = S
P_IR    =    P_IR/180*pi;
P_Vis1D = P_Vis1D/180*pi;
P_Sig1D = P_Sig1D/180*pi;

E = EPolar3(P_Sig1D,P_Vis1D,P_IR); % => [8,1]

E_J_ResLF_Avg = E * J_ResLF_Avg; % => [1,N]

%% Bin signal
AccuGrid = Bin1D(Ex_F1,E_J_ResLF_Avg,FreqRange);

%% Output
OneDSFG.FilesName    = Structure.FilesName;
OneDSFG.SpecType     = 'SFG';
OneDSFG.Response1D   = AccuGrid;
OneDSFG.freq_OneD    = FreqRange;
OneDSFG.H            = H;
OneDSFG.Mu           = Mu;
OneDSFG.Alpha        = Alpha;
OneDSFG.MolFrame     = ResLF;
OneDSFG.R_Avg        = R_Avg;
OneDSFG.LabFrame     = ResLF_Avg;
OneDSFG.Jones        = J;
OneDSFG.JLabFrame    = J_ResLF_Avg;
OneDSFG.E            = E;
OneDSFG.EJLabFrame   = E_J_ResLF_Avg;
