function S_PDB_AmideI = ConstructAmideIPDB(S_PDB,GUI_Inputs)
% The constructiuon function used the SturctureData object from the PDB
% parser (Load_PDB) with the basic info:
% S_PDB.XYZ       
% S_PDB.AtomName  
% S_PDB.FilesName 
% S_PDB.Extra.PDB

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultPhi_D   = 0;
defaultPsi_D   = 0;
defaultTheta_D = 0;
defaultNLFreq     = 1650;
defaultAnharm     = 10;
defaultLFreq      = 1630;
defaultL_Index    = 'None';

addOptional(INPUT,'Phi_D'   ,defaultPhi_D  );
addOptional(INPUT,'Psi_D'   ,defaultPsi_D  );
addOptional(INPUT,'Theta_D' ,defaultTheta_D);
addOptional(INPUT,'NLFreq'   , defaultNLFreq    );
addOptional(INPUT,'Anharm'   , defaultAnharm    );
addOptional(INPUT,'LFreq'    , defaultLFreq     );
addOptional(INPUT,'L_Index'  , defaultL_Index   );

parse(INPUT,GUI_Inputs_C{:});

Phi_D   = INPUT.Results.Phi_D;
Psi_D   = INPUT.Results.Psi_D;
Theta_D = INPUT.Results.Theta_D;
NLFreq    = INPUT.Results.NLFreq;
Anharm    = INPUT.Results.Anharm;
LFreq     = INPUT.Results.LFreq;
L_Index   = INPUT.Results.L_Index;

%% Add AmideI mode info into StructureData class PDB
S_PDB_AmideI = SD_GetAmideI(S_PDB);

%% Rotate the molecule
Phi   = Phi_D/180*pi;
Psi   = Psi_D/180*pi;
Theta = Theta_D/180*pi;
R     = R1_ZYZ_0(Phi,Psi,Theta);
S_PDB_AmideI = SD_Rot(S_PDB_AmideI,R);

% isotop labeling
LocFreq = ones(S_PDB_AmideI.Nmodes,1).*NLFreq;
if ~isempty(L_Index)
    LocFreq(L_Index) = LFreq;  
end
S_PDB_AmideI.LocFreq   = LocFreq;
S_PDB_AmideI.LocAnharm = ones(S_PDB_AmideI.Nmodes,1).*Anharm;

% Add figure drawing function and GUI inputs
S_PDB_AmideI.GUI_Inputs = GUI_Inputs;
S_PDB_AmideI.hPlotFunc  = @PlotXYZfiles_AmideI;

% Calculate One Exciton Hamiltonian
S_PDB_AmideI = SD_1ExH(S_PDB_AmideI); 