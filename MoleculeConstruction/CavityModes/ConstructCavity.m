function S_C = ConstructCavity(GUI_Inputs)
%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultNCModes      = 1;
defaultFSR          = 100;
defaultNLFreq       = 1650;
defaultLocFreqType  = 'Cavity';
defaultCouplingType = 'Zero';

% add Optional inputs / Parameters
addOptional(INPUT,'NCModes'     ,defaultNCModes);
addOptional(INPUT,'FSR'         ,defaultFSR);
addOptional(INPUT,'NLFreq'      ,defaultNLFreq);
addOptional(INPUT,'LocFreqType' ,defaultLocFreqType);
addOptional(INPUT,'CouplingType',defaultCouplingType);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
NCModes      = INPUT.Results.NCModes;
FSR          = INPUT.Results.FSR;
NLFreq       = INPUT.Results.NLFreq;
LocFreqType  = INPUT.Results.LocFreqType;
CouplingType = INPUT.Results.CouplingType;

%% Constrcut local frequencies array
FSR_Array    = (0:FSR:(NCModes-1)*FSR)';
LocFreqArray = ones(NCModes,1)*NLFreq + FSR_Array;

%% Construct essential structure data of the cavity
% mu and alpha of carboxylic acid calculated by G09
% I use these value as a reference to decide a meaningful magnitude 
% mu_Mol    = [-5.43155E+00,  1.653500E+01, -3.17990E-05];
% norm(mu_Mol) = 17.4043 [KM/mol]
% unit conversion: [KM/mol] * 0.23344 = [Debye] 
mu_Mol = [17.4043, 0.0000, 0.0000].*10; % 10 times larger than the carboxylic acid

S_C = StructureData;
S_C.XYZ       = zeros(NCModes,3);
S_C.AtomName  = cellstr(repmat('C',NCModes,1)); % place carbon atoms at the center of the cavity mode for the CoM method to work
S_C.LocCenter = zeros(NCModes,3);
S_C.LocFreq   = LocFreqArray;
S_C.LocMu     = repmat(mu_Mol,NCModes,1);
S_C.LocAlpha  = zeros(NCModes,9); % let's assume ther's no raman active cavity modes

%% Feeds the default GUI_Options back into the StructureData 
GUI_Inputs.CouplingType = CouplingType;
GUI_Inputs.LocFreqType  = LocFreqType;

%% Others
S_C.GUI_Inputs = GUI_Inputs;
S_C.FilesName = [num2str(NCModes),' cavity modes'];
S_C.hPlotFunc = @PlotCavity;
S_C.Extra.CavityInd = 1:NCModes; % add index for comb2 to tell which modes are cavity modes

%% Calculate One Exciton Hamiltonian
S_C = SD_1ExH(S_C); 