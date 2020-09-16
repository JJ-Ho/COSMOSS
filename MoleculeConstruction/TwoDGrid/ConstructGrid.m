function S_V1N_V2M = ConstructGrid(S_Monomer)
% This function form a 2D grid of the input monomer

%% Inputs parser
GUI_Inputs = S_Monomer.GUI_Inputs;
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultDelta_Phi       = 0;
defaultDelta_Psi       = 0;
defaultDelta_Theta     = 0;
defaultVec_1           = [7,0,0];
defaultVec_2           = [0,4,0];
defaultN_1             = 2;
defaultN_2             = 3;

% add Optional inputs / Parameters
addOptional(INPUT,'Delta_Phi'      ,defaultDelta_Phi      );
addOptional(INPUT,'Delta_Psi'      ,defaultDelta_Psi      );
addOptional(INPUT,'Delta_Theta'    ,defaultDelta_Theta    );
addOptional(INPUT,'Vec_1'          ,defaultVec_1          );
addOptional(INPUT,'Vec_2'          ,defaultVec_2          );
addOptional(INPUT,'N_1'            ,defaultN_1            );
addOptional(INPUT,'N_2'            ,defaultN_2            );

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
Vec_1           = INPUT.Results.Vec_1          ;
Vec_2           = INPUT.Results.Vec_2          ;
N_1             = INPUT.Results.N_1            ;
N_2             = INPUT.Results.N_2            ;

%% Set the molecular frame SD to Lab frame
S_Monomer = R_MF2LF(S_Monomer);

%% Extend the monomer along the two directions
S_V1N     = SD_TransN(S_Monomer,Vec_1,N_1);
S_V1N_V2M = SD_TransN(S_V1N,Vec_2,N_2);

%% Deal with files name 
S_V1N_V2M.FilesName = [ '2D_Grid_V1' num2str(N_1) '_V2' num2str(N_2)];

%% Save the Translational copy info into Extra and Children
Extra        = S_V1N_V2M.Extra;
Extra.N_Vec1 = N_1;
Extra.N_Vec2 = N_2;
Extra.Vec_1  = Vec_1;
Extra.Vec_2  = Vec_2;

S_V1N_V2M.Extra    = Extra;

%% Calculate One Exciton Hamiltonian
S_V1N_V2M = SD_1ExH(S_V1N_V2M); 
