function  [SpectraGrid,Response] = TwoDSFG_Main(PDB_Data,GUI_Inputs)
%% TwoDSFG_AmideI
%  
%   Given a initial stucture (pdb), this script will simulate its 2DSFG
%   spectrum.
% 

% ------- Version log -----------------------------------------------------
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
% clear all
% 
% PDB_Data = GetAcid;
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
defaultBeta_NN      = 0.8; % 0.8 cm-1 according to Lauren's PNAS paper (doi/10.1073/pnas.1117704109); that originate from Min Cho's paper (doi:10.1063/1.1997151)
defaultAvg_Phi      = 0;
defaultAvg_Theta    = 0;
defaultAvg_Psi      = 0;
defaultAvg_Rot      = 1;
defaultAvg_Mirror   = 1;
defaultA_Pump       = 90;
defaultA_Probe      = 90;
defaultA_Vis2D      = 90;
defaultA_Sig2D      = 90;
defaultP_Pump1      = 0;
defaultP_Pump2      = 0;
defaultP_Probe      = 0;
defaultP_Vis2D      = 0;
defaultP_Sig2D      = 0;

addOptional(INPUT,'FreqRange'   ,defaultFreqRange);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);
addOptional(INPUT,'Avg_Phi'     ,defaultAvg_Phi);
addOptional(INPUT,'Avg_Theta'   ,defaultAvg_Theta);
addOptional(INPUT,'Avg_Psi'     ,defaultAvg_Psi);
addOptional(INPUT,'Avg_Rot'     ,defaultAvg_Rot);
addOptional(INPUT,'Avg_Mirror'  ,defaultAvg_Mirror);
addOptional(INPUT,'A_Pump'      ,defaultA_Pump);
addOptional(INPUT,'A_Probe'     ,defaultA_Probe);
addOptional(INPUT,'A_Vis2D'     ,defaultA_Vis2D);
addOptional(INPUT,'A_Sig2D'     ,defaultA_Sig2D);
addOptional(INPUT,'P_Pump1'     ,defaultP_Pump1);
addOptional(INPUT,'P_Pump2'     ,defaultP_Pump2);
addOptional(INPUT,'P_Probe'     ,defaultP_Probe);
addOptional(INPUT,'P_Vis2D'     ,defaultP_Vis2D);
addOptional(INPUT,'P_Sig2D'     ,defaultP_Sig2D);
addOptional(INPUT,'CouplingType',defaultCouplingType);


parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
FreqRange    = INPUT.Results.FreqRange;
CouplingType = INPUT.Results.CouplingType;
Beta_NN      = INPUT.Results.Beta_NN;
Avg_Phi      = INPUT.Results.Avg_Phi;
Avg_Theta    = INPUT.Results.Avg_Theta;
Avg_Psi      = INPUT.Results.Avg_Psi;
A_Pump       = INPUT.Results.A_Pump;
A_Probe      = INPUT.Results.A_Probe;
A_Vis2D      = INPUT.Results.A_Vis2D;
A_Sig2D      = INPUT.Results.A_Sig2D;
P_Pump1      = INPUT.Results.P_Pump1;
P_Pump2      = INPUT.Results.P_Pump2;
P_Probe      = INPUT.Results.P_Probe;
P_Vis2D      = INPUT.Results.P_Vis2D;
P_Sig2D      = INPUT.Results.P_Sig2D;
Avg_Rot      = INPUT.Results.Avg_Rot;
Avg_Mirror   = INPUT.Results.Avg_Mirror; 

%% Correct the units

% Orientation = Orientation/180*pi; % turn to radius unit
Avg_Phi_R   =   Avg_Phi/180*pi;
Avg_Psi_R   =   Avg_Psi/180*pi;
Avg_Theta_R = Avg_Theta/180*pi;

% Turn degrees into radius
A_Pump  =  A_Pump /180*pi;
A_Probe =  A_Probe/180*pi;
A_Vis2D =  A_Vis2D/180*pi;
A_Sig2D =  A_Sig2D/180*pi;

% Polarization Angles of incident beams
P_Pump1 = P_Pump1/180*pi;   
P_Pump2 = P_Pump2/180*pi;
P_Probe = P_Probe/180*pi;
P_Vis2D = P_Vis2D/180*pi;
P_Sig2D = P_Sig2D/180*pi;
    
%% Call TwoExcitonH to calculate H,mu and alpha under exciton basis

H = ExcitonH(PDB_Data,...
             'ExMode'  ,'TwoEx',...
             'CouplingType',CouplingType,...
             'Beta_NN' ,Beta_NN);

Mu    = MuAlphaGen(PDB_Data,H,'Mode','Mu');
Alpha = MuAlphaGen(PDB_Data,H,'Mode','Alpha');

Sort_Ex_Freq = H.Sort_Ex_Freq;
Mu_Ex        = Mu.Trans_Ex;
Alpha_Ex     = Alpha.Trans_Ex;

%% Generate Feynman pathway for 2DSFG
Num_Modes = PDB_Data.Num_Modes;
Response  = Feynman_2DSFG_kron(Num_Modes,Sort_Ex_Freq,Alpha_Ex,Mu_Ex);

Response.H = H;
Response.Mu = Mu;
Response.Alpha = Alpha;
%% Decide what kinds of rod rotation average is

switch Avg_Rot
    
    case 1 %'Phi' C_Inf

        R_Avg = R5_ZYZ_1(Avg_Psi_R,Avg_Theta_R);
        
    case 2 %'Psi'

        R_Avg = R5_ZYZ_2(Avg_Phi_R,Avg_Theta_R);
        
    case 3 %'{Phi,Psi}'

        R_Avg = R5_ZYZ_12(Avg_Theta_R);
        
    case 4 %'Isotropic'
        
        R_Avg = R5_ZYZ_123;
    
    case 5 %'No Average'
       
        R_Avg = R5_ZYZ_0(Avg_Phi_R,Avg_Psi_R,Avg_Theta_R);
        
    otherwise
        disp('Avg_Option is not support, dont know how to apply Rotational average...')
end

% Decide Mirror planes

switch Avg_Mirror
    
    case 1 % no mirror plane
        V = [1;1;1];
        
        Mirror_Mask = kron(kron(kron(kron(V,V),V),V),V);
        
    case 2 % sigma v, X=-X, Y=-Y
        V1 = [-1; 1;1];
        V2 = [ 1;-1;1];
        
        Sigma_X = kron(kron(kron(kron(V1,V1),V1),V1),V1);
        Sigma_Y = kron(kron(kron(kron(V2,V2),V2),V2),V2);
        Sigma_X(eq(Sigma_X,-1)) = 0;
        Sigma_Y(eq(Sigma_Y,-1)) = 0;
        
        Mirror_Mask = and(Sigma_X,Sigma_Y);
end

%% Applied rotational avergae on Response in molecular frame

AR1  = zeros(size(Response.R1));
AR2  = zeros(size(Response.R2));
AR3  = zeros(size(Response.R3));
NAR1 = zeros(size(Response.NR1));
NAR2 = zeros(size(Response.NR2));
NAR3 = zeros(size(Response.NR3));

AR1(:,1:3)  = Response.R1(:,1:3);
AR2(:,1:3)  = Response.R2(:,1:3);
AR3(:,1:3)  = Response.R3(:,1:3);
NAR1(:,1:3) = Response.NR1(:,1:3);
NAR2(:,1:3) = Response.NR2(:,1:3);
NAR3(:,1:3) = Response.NR3(:,1:3);

AR1(:,4:end)  = (bsxfun(@times,R_Avg*Response.R1(:,4:end)' ,Mirror_Mask))';
AR2(:,4:end)  = (bsxfun(@times,R_Avg*Response.R2(:,4:end)' ,Mirror_Mask))';
AR3(:,4:end)  = (bsxfun(@times,R_Avg*Response.R3(:,4:end)' ,Mirror_Mask))';
NAR1(:,4:end) = (bsxfun(@times,R_Avg*Response.NR1(:,4:end)',Mirror_Mask))';
NAR2(:,4:end) = (bsxfun(@times,R_Avg*Response.NR2(:,4:end)',Mirror_Mask))';
NAR3(:,4:end) = (bsxfun(@times,R_Avg*Response.NR3(:,4:end)',Mirror_Mask))';

Response.AR1 = AR1;
Response.AR2 = AR2;
Response.AR3 = AR3;
Response.NAR1 = NAR1;
Response.NAR2 = NAR2;
Response.NAR3 = NAR3;

%% Jones Matrix convert XYZ to PS frame

% Note: When I generate J, I aasume the two pump beam has the same input angle
%       That's why here only has one pump input 
J = JonesRef5(A_Sig2D,A_Vis2D,A_Probe,A_Pump);

% 35 =3 freq + 2^5 (ppppp - sssss) of signal
JAR1  = zeros(size(Response.AR1,1),35);
JAR2  = zeros(size(Response.AR2,1),35);
JAR3  = zeros(size(Response.AR3,1),35);
JNAR1 = zeros(size(Response.NAR1,1),35);
JNAR2 = zeros(size(Response.NAR2,1),35);
JNAR3 = zeros(size(Response.NAR3,1),35);

JAR1(:,1:3)  = Response.AR1(:,1:3);
JAR2(:,1:3)  = Response.AR2(:,1:3);
JAR3(:,1:3)  = Response.AR3(:,1:3);
JNAR1(:,1:3) = Response.NAR1(:,1:3);
JNAR2(:,1:3) = Response.NAR2(:,1:3);
JNAR3(:,1:3) = Response.NAR3(:,1:3);

JAR1(:,4:end)  = (J*Response.AR1(:,4:end)')';
JAR2(:,4:end)  = (J*Response.AR2(:,4:end)')';
JAR3(:,4:end)  = (J*Response.AR3(:,4:end)')';
JNAR1(:,4:end) = (J*Response.NAR1(:,4:end)')';
JNAR2(:,4:end) = (J*Response.NAR2(:,4:end)')';
JNAR3(:,4:end) = (J*Response.NAR3(:,4:end)')';

Response.JAR1 = JAR1;
Response.JAR2 = JAR2;
Response.JAR3 = JAR3;
Response.JNAR1 = JNAR1;
Response.JNAR2 = JNAR2;
Response.JNAR3 = JNAR3;

%% E part, Plarization of each incident beams

E = EPolar5(P_Sig2D,P_Vis2D,P_Probe,P_Pump2,P_Pump1);

% 4 = 3 freq + 1 signal
EJAR1  = zeros(size(Response.JAR1,1),4);
EJAR2  = zeros(size(Response.JAR2,1),4);
EJAR3  = zeros(size(Response.JAR3,1),4);
EJNAR1 = zeros(size(Response.JNAR1,1),4);
EJNAR2 = zeros(size(Response.JNAR2,1),4);
EJNAR3 = zeros(size(Response.JNAR3,1),4);

EJAR1(:,1:3)  = Response.JAR1(:,1:3);
EJAR2(:,1:3)  = Response.JAR2(:,1:3);
EJAR3(:,1:3)  = Response.JAR3(:,1:3);
EJNAR1(:,1:3) = Response.JNAR1(:,1:3);
EJNAR2(:,1:3) = Response.JNAR2(:,1:3);
EJNAR3(:,1:3) = Response.JNAR3(:,1:3);

EJAR1(:,4:end)  = (E*Response.JAR1(:,4:end)')';
EJAR2(:,4:end)  = (E*Response.JAR2(:,4:end)')';
EJAR3(:,4:end)  = (E*Response.JAR3(:,4:end)')';
EJNAR1(:,4:end) = (E*Response.JNAR1(:,4:end)')';
EJNAR2(:,4:end) = (E*Response.JNAR2(:,4:end)')';
EJNAR3(:,4:end) = (E*Response.JNAR3(:,4:end)')';

Response.EJAR1 = EJAR1;
Response.EJAR2 = EJAR2;
Response.EJAR3 = EJAR3;
Response.EJNAR1 = EJNAR1;
Response.EJNAR2 = EJNAR2;
Response.EJNAR3 = EJNAR3;

%% Binning of spectra
Response.BinR1  = EJAR1;
Response.BinR2  = EJAR2;
Response.BinR3  = EJAR3;
Response.BinNR1 = EJNAR1;
Response.BinNR2 = EJNAR2;
Response.BinNR3 = EJNAR3;

SpectraGrid = Bin2D(Response,FreqRange);

