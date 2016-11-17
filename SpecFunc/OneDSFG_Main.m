function  OneDSFG = OneDSFG_Main(Structure_Data,GUI_Inputs)
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
defaultFreqRange    = 1650:1750;

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
addOptional(INPUT,'FreqRange'   ,defaultFreqRange);

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
FreqRange    = INPUT.Results.FreqRange;

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

% vectorized version
Alpha_Ex = squeeze(Alpha_Ex(2:Num_Modes+1,1,:))'; %=> [9*N]
Mu_Ex    = squeeze(   Mu_Ex(1,2:Num_Modes+1,:))'; %=> [3*N]
[Mu_Ind,Alpha_Ind] = ndgrid(1:3,1:9);

ResLF = Alpha_Ex(Alpha_Ind(:),:).* Mu_Ex(Mu_Ind(:),:); %=> [27*N]
Ex_Freq = Ex_Freq(2:Num_Modes+1); % => [N*1]

%% Decide what kinds of ensemble average

Dimension = 3; % for SFG

switch Avg_Rot
        
    case 1 %'Phi' C_Inf

        R_Avg = LabFrameAvg('C4',Dimension);
        
    
    case 5 %'No Average'
       
        R_Avg = LabFrameAvg('C1',Dimension);
    
        
    otherwise
        disp('Avg_Angle is not support, dont know how to apply Rotational average...')
        return
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

ResLF_Avg = bsxfun(@times,R_Avg*ResLF,Mirror_Mask); %=> [27*N]

%% Jones Matrix convert XYZ to PS frame
% 
% JLabFrame = [freq, ppp, pps, psp, pss, spp, sps, ssp, sss]
% 
% Turn degrees into radius
A_IR  =  A_IR/180*pi;
A_Vis = A_Vis/180*pi;
A_Sum = A_Sum/180*pi;

J = JonesRef3(A_Sum,A_Vis,A_IR); % => [27,8]

J_ResLF_Avg = J * ResLF_Avg; % => [8,N]

%% E part, Plarization of each incident beams

E = EPolar3(P_Sum,P_Vis,P_IR); % => [8,1]

E_J_ResLF_Avg = E * J_ResLF_Avg; % => [1,N]

%% Bin signal
AccuGrid = Bin1D(Ex_Freq,E_J_ResLF_Avg,FreqRange);

%% Output
OneDSFG.H            = H;
OneDSFG.Mu           = Mu;
OneDSFG.Alpha        = Alpha;
OneDSFG.Num_Modes    = Num_Modes;
OneDSFG.MolFrame     = ResLF;
OneDSFG.R_Avg        = R_Avg;
OneDSFG.LabFrame     = ResLF_Avg;
OneDSFG.Jones        = J;
OneDSFG.JLabFrame    = J_ResLF_Avg;
OneDSFG.E            = E;
OneDSFG.EJLabFrame   = E_J_ResLF_Avg;
OneDSFG.FilesName    = Structure_Data.FilesName;
OneDSFG.CouplingType = CouplingType;
OneDSFG.SpecType     = 'SFG';
OneDSFG.Response1D   = AccuGrid;
OneDSFG.freq_OneD    = FreqRange;
