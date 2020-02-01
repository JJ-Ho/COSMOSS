function Output= ExcitonH(SData,GUI_Inputs)
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
defaultBeta_NN      = 0.8; % 0.8 cm-1 according to Lauren's PNAS paper (doi/10.1073/pnas.1117704109); that originate from Min Cho's paper (doi:10.1063/1.1997151)

addOptional(INPUT,'LocFreqType' ,defaultLocFreqType);
addOptional(INPUT,'CouplingType',defaultCouplingType);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);

parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
LocFreqType  = INPUT.Results.LocFreqType;
CouplingType = INPUT.Results.CouplingType;
Beta_NN      = INPUT.Results.Beta_NN;

%% para variables
% Reassign variable names from StrucInfo
Nmodes  = SData.Nmodes;
LocFreq = SData.LocFreq;

%% Apply Jansen map if needed
if eq(LocFreqType,2)
    [~,dF_Jansen] = Coupling_Jansen(SData);
    LocFreq = LocFreq + dF_Jansen;
end

Beta = Coupling(SData,CouplingType,Beta_NN); % Coupling

%% One Exciton full Hamiltonian
% Zero exciton part of full Hamiltonain 
ZeroExPart = 0;

% The result is in cm-1 unit
OneExPart = bsxfun(@times,eye(Nmodes),LocFreq) + Beta;

% Construct Full H
OneExH = blkdiag(ZeroExPart,OneExPart);

%% Output Variables
Output.Beta   = Beta;
Output.OneExH = OneExH;