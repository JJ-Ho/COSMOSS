function  [SpectraGrid,Response] = TwoDSFG_Main_Sparse(Structure,GUI_Inputs)
%% TwoDSFG_AmideI
%  
%   Given a initial stucture (pdb), this script will simulate its 2DSFG
%   spectrum.
% 

% ------- Version log -----------------------------------------------------
% 
% Ver. 2.0  161127  Reduce and break the ouput of H, Mu, Alpha into 0->1
%                   and 1->2.. This simplify the indexing of matrix element
%                   and helps in future cut-off speed-up upgrade. 
% 
% Ver. 1.3  141016  Add "PolarAng", "INCAng", "RotationalAvg", "Mirror_Plane"
%                   into optional input. Take out unecessary input "handles"
% 
% Ver. 1.2  140717  Add frequency range readin function
% 
% Ver. 1.1  140608  Functionalize H, Mu Alpha
% 
% Ver. 1.0  140607  Copy from previous devolope version  
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
defaultAvg_Rot      = 1;
defaultAvg_Mirror   = 1;
defaultA_Pump1      = 90;
defaultA_Pump2      = 90;
defaultA_Probe      = 90;
defaultA_Vis2D      = 90;
defaultA_Sig2D      = 90;
defaultP_Pump1      = 0;
defaultP_Pump2      = 0;
defaultP_Probe      = 0;
defaultP_Vis2D      = 0;
defaultP_Sig2D      = 0;
defaultPCutOff      = 0;

addOptional(INPUT,'FreqRange'   ,defaultFreqRange);
addOptional(INPUT,'Avg_Rot'     ,defaultAvg_Rot);
addOptional(INPUT,'Avg_Mirror'  ,defaultAvg_Mirror);
addOptional(INPUT,'A_Pump1'     ,defaultA_Pump1);
addOptional(INPUT,'A_Pump2'     ,defaultA_Pump2);
addOptional(INPUT,'A_Probe'     ,defaultA_Probe);
addOptional(INPUT,'A_Vis2D'     ,defaultA_Vis2D);
addOptional(INPUT,'A_Sig2D'     ,defaultA_Sig2D);
addOptional(INPUT,'P_Pump1'     ,defaultP_Pump1);
addOptional(INPUT,'P_Pump2'     ,defaultP_Pump2);
addOptional(INPUT,'P_Probe'     ,defaultP_Probe);
addOptional(INPUT,'P_Vis2D'     ,defaultP_Vis2D);
addOptional(INPUT,'P_Sig2D'     ,defaultP_Sig2D);
addOptional(INPUT,'PCutOff'     ,defaultPCutOff);

parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
FreqRange    = INPUT.Results.FreqRange;
Avg_Rot      = INPUT.Results.Avg_Rot;
Avg_Mirror   = INPUT.Results.Avg_Mirror; 
A_Pump1      = INPUT.Results.A_Pump1;
A_Pump2      = INPUT.Results.A_Pump2;
A_Probe      = INPUT.Results.A_Probe;
A_Vis2D      = INPUT.Results.A_Vis2D;
A_Sig2D      = INPUT.Results.A_Sig2D;
P_Pump1      = INPUT.Results.P_Pump1;
P_Pump2      = INPUT.Results.P_Pump2;
P_Probe      = INPUT.Results.P_Probe;
P_Vis2D      = INPUT.Results.P_Vis2D;
P_Sig2D      = INPUT.Results.P_Sig2D;
PCutOff      = INPUT.Results.PCutOff;

%% Call TwoExcitonH to calculate H,mu and alpha under exciton basis
H = ExcitonH(Structure,GUI_Inputs,'TwoEx');    

Mu    = MuAlphaGen(Structure,H,'Mode','Mu');
Alpha = MuAlphaGen(Structure,H,'Mode','Alpha');

Ex_F1   = H.Sort_Ex_F1;
Ex_F2   = H.Sort_Ex_F2;
M_Ex_01 = Mu.M_Ex_01;
M_Ex_12 = Mu.M_Ex_12;
A_Ex_01 = Alpha.M_Ex_01;
A_Ex_12 = Alpha.M_Ex_12;

%% Decide what kinds of rod rotation average is
% note:sparse version dose not have mirror plane impemented, yet...
[R_Avg,Mirror_Mask,~,~] = LabFrameAvg(Avg_Rot,Avg_Mirror,5);


%% Jones Matrix convert XYZ to PS frame
% Laser incident angles between laser beam and surface normal.
A_Pump1 = A_Pump1/180*pi;
A_Pump2 = A_Pump2/180*pi;
A_Probe = A_Probe/180*pi;
A_Vis2D = A_Vis2D/180*pi;
A_Sig2D = A_Sig2D/180*pi;

J = JonesRef5(A_Sig2D,A_Vis2D,A_Probe,A_Pump2,A_Pump1); % Take [radius]

%% E part, Plarization of each incident beams
% Polarization Angles of incident beams, 0 = P, 90 = S
P_Pump1 = P_Pump1/180*pi;   
P_Pump2 = P_Pump2/180*pi;
P_Probe = P_Probe/180*pi;
P_Vis2D = P_Vis2D/180*pi;
P_Sig2D = P_Sig2D/180*pi;

E = EPolar5(P_Sig2D,P_Vis2D,P_Probe,P_Pump2,P_Pump1); % Take [radius]

%% Generate Feynman pathway for 2DSFG
EJR = E*J*R_Avg;

[SpectraGrid,Freq,Int,Index,CutOff] = Feynman_2DSFG_Vec_Sparse(PCutOff,...
                                                        FreqRange,...
                                                        EJR,...
                                                        Ex_F1,...
                                                        Ex_F2,...
                                                        A_Ex_01,...
                                                        A_Ex_12,...
                                                        M_Ex_01,...
                                                        M_Ex_12);

%% Group up other outputs
Response.H = H;
Response.Mu = Mu;
Response.Alpha = Alpha;

Response.Freq   = Freq;
Response.Int    = Int;
Response.Index  = Index;
Response.CutOff = CutOff;

Response.EJR = EJR;
Response.SpecType = '2DSFG';
