function Output= ExcitonH(Structure,GUI_Inputs)
%% TwoExcitonH
%  
% Given StrucInfo that generate by GetAmindeI, this script generate One or
% Two Exciton Hamiltonian and also export eigenvalues and eigenvectors for
% generating transition dipole and Raman tensor matrix. 
% 
% Todo: fix [Improve] part
%       Integrate Diagonal disorder part in.
% 
% Copyright Jia-Jung Ho, 2013-2020

%% Debug
% Structure  = Data_COSMOSS.Structure;
% GUI_Inputs = ParseGUI_Main(Data_COSMOSS.hGUIs);

%% Inputs parser
% Turn Output from Read GUI to cell array
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

defaultLocFreqType  = 1;
defaultCouplingType = 'TDC';
defaultSampling     = 0;
defaultP_FlucCorr   = 100;
% defaultDD_FWHM      = 0;
defaultODD_FWHM     = 0;
defaultBeta_NN      = 0.8; % 0.8 cm-1 according to Lauren's PNAS paper (doi/10.1073/pnas.1117704109); that originate from Min Cho's paper (doi:10.1063/1.1997151)

addOptional(INPUT,'LocFreqType' ,defaultLocFreqType);
addOptional(INPUT,'CouplingType',defaultCouplingType);
addOptional(INPUT,'Sampling'    ,defaultSampling);
addOptional(INPUT,'P_FlucCorr'  ,defaultP_FlucCorr);
% addOptional(INPUT,'DD_FWHM'     ,defaultDD_FWHM);
addOptional(INPUT,'ODD_FWHM'    ,defaultODD_FWHM);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);

parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
LocFreqType  = INPUT.Results.LocFreqType;
CouplingType = INPUT.Results.CouplingType;
Sampling     = INPUT.Results.Sampling;
P_FlucCorr   = INPUT.Results.P_FlucCorr;
% DD_FWHM      = INPUT.Results.DD_FWHM;
ODD_FWHM     = INPUT.Results.ODD_FWHM;
Beta_NN      = INPUT.Results.Beta_NN;

%% para variables
% Reassign variable names from StrucInfo
Nmodes          = Structure.Nmodes;
LocFreq         = Structure.LocFreq;
DiagDisorder    = Structure.DiagDisorder;
% OffDiagDisorder = Structure.OffDiagDisorder; have not implement yet

% check if apply random sampling
if Sampling
    DD_std  = DiagDisorder./(2*sqrt(2*log(2)));
    ODD_std = ODD_FWHM/(2*sqrt(2*log(2)));
else
    DD_std  = 0;
    ODD_std = 0;
end

%% Diagonal disorder if any
if eq(LocFreqType,2)
    [~,dF_Jansen] = Coupling_Jansen(Structure);
    LocFreq = LocFreq + dF_Jansen;
end

P_FlucCorr = P_FlucCorr/100; % turn percentage to number within 0~1

Correlation_Dice = rand;
if Correlation_Dice < P_FlucCorr
    dF_DD = DD_std.*(randn(1,1).*ones(Nmodes,1));
else 
    dF_DD = DD_std.*randn(Nmodes,1); 
end
LocFreq = LocFreq + dF_DD;

%% Off diagonal disorder
dBeta   = ODD_std*randn(Nmodes);
dBeta   = (dBeta + dBeta')./2; % symetrize
Beta    = Coupling(Structure,CouplingType,Beta_NN); % Coupling
Beta    = Beta + dBeta;

%% Zero exciton part of full Hamiltonain 
ZeroExPart = 0;

%% One Exciton part of full Hamiltonian
% The result is in cm-1 unit
OneExPart = bsxfun(@times,eye(Nmodes),LocFreq) + Beta;

% Construct Full H
H = blkdiag(ZeroExPart,OneExPart);

%% Output Variables
Output.Nmodes        = Nmodes;
Output.Beta          = Beta;
Output.H             = H;
Output.OneExPart     = OneExPart;