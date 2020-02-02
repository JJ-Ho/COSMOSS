function obj_1ExH = SD_1ExH(obj_SD)
%% OneExcitonH generation method of StructureData class
%  
% Given StructureData, this script generate 1 Exciton Hamiltonian
% Copyright Jia-Jung Ho, 2013-2020

%% GUI Inputs parser
GUI_Inputs = obj_SD.GUI_Inputs;

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

%% Apply Jansen map if needed
LocFreq = obj_SD.LocFreq;
if strcmp(LocFreqType,'Jensen Map')
    [~,dF_Jansen] = Coupling_Jansen(obj_SD);
    LocFreq = LocFreq + dF_Jansen;
end

Beta = Coupling(obj_SD,CouplingType,Beta_NN); % Coupling

%% Output Variables
% Note: the One exciton Hamiltonian will be calculated providing the
% following data exisited in the StructureData
obj_1ExH         = obj_SD;
obj_1ExH.LocFreq = LocFreq;
obj_1ExH.Beta    = Beta;