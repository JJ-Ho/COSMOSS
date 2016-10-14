function  OneDSFG = OneDSFG_Main(Structure_Data,GUI_Inputs)
%% TwoDSFG_AmideI_Main

% ------- Version log -----------------------------------------------------
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
defaultAvg_Phi      = 0;
defaultAvg_Theta    = 0;
defaultAvg_Psi      = 0;
defaultAvg_Rot      = 1;
defaultAvg_Mirror   = 1;
defaultA_IR         = 90;
defaultA_Vis        = 90;
defaultA_Sum        = 90;
defaultP_IR         = 90;
defaultP_Vis        = 90;
defaultP_Sum        = 90;
defaultBeta_NN      = 0.8;

% add Optional inputs / Parameters
addOptional(INPUT,'CouplingType',defaultCouplingType);
addOptional(INPUT,'Avg_Phi'     ,defaultAvg_Phi);
addOptional(INPUT,'Avg_Theta'   ,defaultAvg_Theta);
addOptional(INPUT,'Avg_Psi'     ,defaultAvg_Psi);
addOptional(INPUT,'Avg_Rot'     ,defaultAvg_Rot);
addOptional(INPUT,'Avg_Mirror'  ,defaultAvg_Mirror);
addOptional(INPUT,'A_IR'        ,defaultA_IR);
addOptional(INPUT,'A_Vis'       ,defaultA_Vis);
addOptional(INPUT,'A_Sum'       ,defaultP_Sum);
addOptional(INPUT,'P_IR'        ,defaultP_IR);
addOptional(INPUT,'P_Vis'       ,defaultP_Vis);
addOptional(INPUT,'P_Sum'       ,defaultA_Sum);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names

CouplingType = INPUT.Results.CouplingType;
Avg_Phi      = INPUT.Results.Avg_Phi;
Avg_Theta    = INPUT.Results.Avg_Theta;
Avg_Psi      = INPUT.Results.Avg_Psi;
Avg_Rot      = INPUT.Results.Avg_Rot;
Avg_Mirror   = INPUT.Results.Avg_Mirror;
A_IR         = INPUT.Results.A_IR;
A_Vis        = INPUT.Results.A_Vis;
A_Sum        = INPUT.Results.A_Sum;
P_IR         = INPUT.Results.P_IR;
P_Vis        = INPUT.Results.P_Vis;
P_Sum        = INPUT.Results.P_Sum;
Beta_NN      = INPUT.Results.Beta_NN;

%% Call OneExcitonH to calculate H,mu and alpha under exciton basis

H = ExcitonH(Structure_Data,'ExMode','OneEx','CouplingType',CouplingType,'Beta_NN',Beta_NN);

Ex_Freq = H.Sort_Ex_Freq;

% construct mu,alpha
Mu    = MuAlphaGen(Structure_Data,H,'Mode','Mu');
Alpha = MuAlphaGen(Structure_Data,H,'Mode','Alpha');

Mu_Ex    = Mu.Trans_Ex;
Alpha_Ex = Alpha.Trans_Ex;

%% Generate Molecular frame SFG Responses
Num_Modes = Structure_Data.Num_Modes;

ResMolFrame = zeros(Num_Modes,3^3+1);

for N = 1:Num_Modes
    ResMolFrame(N,1)     = Ex_Freq(N+1);
    ResMolFrame(N,2:end) = kron(squeeze(Alpha_Ex(N+1,1,:)),squeeze(Mu_Ex(1,N+1,:)));
end

%% Decide what kinds of ensemble average

% Orientation = Orientation/180*pi; % turn to radius unit
Avg_Phi_R   =   Avg_Phi/180*pi;
Avg_Psi_R   =   Avg_Psi/180*pi;
Avg_Theta_R = Avg_Theta/180*pi;

switch Avg_Rot
        
    case 1 %'Phi' C_Inf

        R_Avg = R3_ZYZ_1(Avg_Psi_R,Avg_Theta_R);
        
    case 2 %'Psi'

        R_Avg = R3_ZYZ_2(Avg_Phi_R,Avg_Theta_R);
        
    case 3 %'{Phi,Psi}'

        R_Avg = R3_ZYZ_12(Avg_Theta_R);
        
    case 4 %'Isotropic'
        
        R_Avg = R3_ZYZ_123;
    
    case 5 %'No Average'
       
        R_Avg = R3_ZYZ_0(Avg_Phi_R,Avg_Psi_R,Avg_Theta_R);
    
        
    otherwise
        disp('Avg_Angle is not support, dont know how to apply Rotational average...')
end

% Decide Mirror planes
switch Avg_Mirror
    
    case 1 % no mirror plane
        V = [1;1;1];
        
        Mirror_Mask = kron(kron(V,V),V);
        
    case 2 % sigma v, X=-X, Y=-Y
        V1 = [-1; 1;1];
        V2 = [ 1;-1;1];
        
        Sigma_X = kron(kron(V1,V1),V1);
        Sigma_Y = kron(kron(V2,V2),V2);
        Sigma_X(eq(Sigma_X,-1)) = 0;
        Sigma_Y(eq(Sigma_Y,-1)) = 0;
        
        Mirror_Mask = and(Sigma_X,Sigma_Y);
end

%% Applied ensemble avergae on Responses in molecular frame

ResLabFrame  = zeros(size(ResMolFrame));

ResLabFrame(:,1) = ResMolFrame(:,1); 
ResLabFrame(:,2:end) = (bsxfun(@times,R_Avg*ResMolFrame(:,2:end)',Mirror_Mask))';

%% Jones Matrix convert XYZ to PS frame
% 
% JLabFrame = [freq, ppp, pps, psp, pss, spp, sps, ssp, sss]
% 
% Turn degrees into radius
A_IR  =  A_IR/180*pi;
A_Vis = A_Vis/180*pi;
A_Sum = A_Sum/180*pi;

J = JonesRef3(A_Sum,A_Vis,A_IR);

JLabFrame  = zeros(size(ResLabFrame,1),9); % 9 = 1 freq + 2^3 (ppp - sss) of signal

JLabFrame(:,1)      = ResLabFrame(:,1);
JLabFrame(:,2:end)  = (J*ResLabFrame(:,2:end)')';

%% E part, Plarization of each incident beams

E = EPolar3(P_Sum,P_Vis,P_IR);

EJLabFrame  = zeros(size(JLabFrame,1),1); % Only one signal compose of all polarization combination

EJLabFrame(:,1) = JLabFrame(:,1);
EJLabFrame(:,2) = (E*JLabFrame(:,2:end)')';

%% Output
OneDSFG.H            = H;
OneDSFG.Mu           = Mu;
OneDSFG.Alpha        = Alpha;
OneDSFG.Num_Modes    = Num_Modes;
OneDSFG.MolFrame     = ResMolFrame;
OneDSFG.R_Avg        = R_Avg;
OneDSFG.LabFrame     = ResLabFrame;
OneDSFG.Jones        = J;
OneDSFG.JLabFrame    = JLabFrame;
OneDSFG.E            = E;
OneDSFG.EJLabFrame   = EJLabFrame;
OneDSFG.FilesName    = Structure_Data.FilesName;
OneDSFG.CouplingType = CouplingType;
OneDSFG.SpecType     = 'SFG';
