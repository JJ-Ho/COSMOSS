function S_C = ConstructCavity(GUI_Inputs)
%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultNCModes  = 1;

% add Optional inputs / Parameters
addOptional(INPUT,'NCModes',defaultNCModes);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
NCModes = INPUT.Results.NCModes;

%% Construct essential structure data of the cavity
S_C = StructureData;
S_C.XYZ       = zeros(NCModes,3);
S_C.AtomName  = cellstr(repmat('C',NCModes,1)); % place carbon atoms at the center of the cavity mode for the CoM method to work
S_C.LocCenter = zeros(NCModes,3);
S_C.LocMu     = zeros(NCModes,3);
S_C.LocAlpha  = zeros(NCModes,9); % raman tensor vector form [N x 9]

%% Decide coupling types 
GUI_Inputs.CouplingType = 'Zero';

%% Others
S_C.GUI_Inputs = GUI_Inputs;
S_C.FilesName = [num2str(NCModes),' cavity modes'];
% hPlotFunc
S_C.Extra.CavityInd = 1:NCModes; % add index for comb2 to tell which modes are cavity modes

%% Calculate One Exciton Hamiltonian
S_C = SD_1ExH(S_C); 