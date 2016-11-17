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

%% Correct the units

% Turn degrees into radius (not work for BoxCard geometry yet)
A_Pump1  =  A_Pump1 /180*pi;
A_Pump2  =  A_Pump2 /180*pi;
A_Probe  =  A_Probe/180*pi;
A_Sig2D  =  A_Sig2D/180*pi;

% Polarization Angles of incident beams
P_Pump1  = P_Pump1/180*pi;
P_Pump2  = P_Pump2/180*pi;
P_Probe  = P_Probe/180*pi;
P_Sig2D  = P_Sig2D/180*pi;
    
%% Call TwoExcitonH to calculate H,mu and alpha under exciton basis
H = ExcitonH(PDB_Data,...
             'ExMode'  ,'TwoEx',...
             'CouplingType',CouplingType,...
             'Beta_NN' ,Beta_NN);

Mu = MuAlphaGen(PDB_Data,H,'Mode','Mu');

Sort_Ex_Freq = H.Sort_Ex_Freq;
Mu_Ex        = Mu.Trans_Ex;

%% Generate Feynman pathway for 2DSFG
Num_Modes = PDB_Data.Num_Modes;
% Response  = Feynman_2DIR_kron(Num_Modes,Sort_Ex_Freq,Mu_Ex);
Response  = Feynman_2DIR_Vec(Num_Modes,Sort_Ex_Freq,Mu_Ex);
Response.H  = H;
Response.Mu = Mu;

%% Decide what kinds of rod rotation average is

R_Avg = LabFrameAvg('Isotropic',4);

%% Applied rotational average on Response in molecular frame

% AR1  = zeros(size(Response.R1));
% AR2  = zeros(size(Response.R2));
% AR3  = zeros(size(Response.R3));
% NAR1 = zeros(size(Response.NR1));
% NAR2 = zeros(size(Response.NR2));
% NAR3 = zeros(size(Response.NR3));
% 
% AR1(:,1:3)  = Response.R1(:,1:3);
% AR2(:,1:3)  = Response.R2(:,1:3);
% AR3(:,1:3)  = Response.R3(:,1:3);
% NAR1(:,1:3) = Response.NR1(:,1:3);
% NAR2(:,1:3) = Response.NR2(:,1:3);
% NAR3(:,1:3) = Response.NR3(:,1:3);

AR1  = R_Avg*Response.R1 ;
AR2  = R_Avg*Response.R2 ;
AR3  = R_Avg*Response.R3 ;
NAR1 = R_Avg*Response.NR1;
NAR2 = R_Avg*Response.NR2;
NAR3 = R_Avg*Response.NR3;

Response.AR1 = AR1;
Response.AR2 = AR2;
Response.AR3 = AR3;
Response.NAR1 = NAR1;
Response.NAR2 = NAR2;
Response.NAR3 = NAR3;

%% Jones Matrix convert XYZ to PS frame

% Note: When I generate J, I aasume the two pump beam has the same input angle
%       That's why here only has one pump input 
J = JonesTrans4(A_Sig2D,A_Probe,A_Pump2,A_Pump1);

% % 19 =3 freq + 2^4 (pppp - ssss) of signal
% JAR1  = zeros(size(Response.AR1,1),19);
% JAR2  = zeros(size(Response.AR2,1),19);
% JAR3  = zeros(size(Response.AR3,1),19);
% JNAR1 = zeros(size(Response.NAR1,1),19);
% JNAR2 = zeros(size(Response.NAR2,1),19);
% JNAR3 = zeros(size(Response.NAR3,1),19);
% 
% JAR1(:,1:3)  = Response.AR1(:,1:3);
% JAR2(:,1:3)  = Response.AR2(:,1:3);
% JAR3(:,1:3)  = Response.AR3(:,1:3);
% JNAR1(:,1:3) = Response.NAR1(:,1:3);
% JNAR2(:,1:3) = Response.NAR2(:,1:3);
% JNAR3(:,1:3) = Response.NAR3(:,1:3);

JAR1  = J*Response.AR1;
JAR2  = J*Response.AR2;
JAR3  = J*Response.AR3;
JNAR1 = J*Response.NAR1;
JNAR2 = J*Response.NAR2;
JNAR3 = J*Response.NAR3;

Response.JAR1 = JAR1;
Response.JAR2 = JAR2;
Response.JAR3 = JAR3;
Response.JNAR1 = JNAR1;
Response.JNAR2 = JNAR2;
Response.JNAR3 = JNAR3;

%% E part, Plarization of each incident beams

E = EPolar4(P_Sig2D,P_Probe,P_Pump2,P_Pump1);

% % 4 = 3 freq + 1 signal
% EJAR1  = zeros(size(Response.JAR1,1),4);
% EJAR2  = zeros(size(Response.JAR2,1),4);
% EJAR3  = zeros(size(Response.JAR3,1),4);
% EJNAR1 = zeros(size(Response.JNAR1,1),4);
% EJNAR2 = zeros(size(Response.JNAR2,1),4);
% EJNAR3 = zeros(size(Response.JNAR3,1),4);
% 
% EJAR1(:,1:3)  = Response.JAR1(:,1:3);
% EJAR2(:,1:3)  = Response.JAR2(:,1:3);
% EJAR3(:,1:3)  = Response.JAR3(:,1:3);
% EJNAR1(:,1:3) = Response.JNAR1(:,1:3);
% EJNAR2(:,1:3) = Response.JNAR2(:,1:3);
% EJNAR3(:,1:3) = Response.JNAR3(:,1:3);

EJAR1  = E*Response.JAR1;
EJAR2  = E*Response.JAR2;
EJAR3  = E*Response.JAR3;
EJNAR1 = E*Response.JNAR1;
EJNAR2 = E*Response.JNAR2;
EJNAR3 = E*Response.JNAR3;

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

