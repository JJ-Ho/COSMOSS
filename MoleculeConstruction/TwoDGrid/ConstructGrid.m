function S_V1N_V2M = ConstructGrid(S_Monomer,GUI_Inputs)
%% Inputs parser
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
defaultNLFreq          = 1720;
defaultAnharm          = 20;
defaultLFreq           = 1700;
defaultL_Index         = 'None';
defaultDiagDisorder    = 0;
defaultOffDiagDisorder = 0;

% add Optional inputs / Parameters
addOptional(INPUT,'Delta_Phi'      ,defaultDelta_Phi      );
addOptional(INPUT,'Delta_Psi'      ,defaultDelta_Psi      );
addOptional(INPUT,'Delta_Theta'    ,defaultDelta_Theta    );
addOptional(INPUT,'Vec_1'          ,defaultVec_1          );
addOptional(INPUT,'Vec_2'          ,defaultVec_2          );
addOptional(INPUT,'N_1'            ,defaultN_1            );
addOptional(INPUT,'N_2'            ,defaultN_2            );
addOptional(INPUT,'NLFreq'         ,defaultNLFreq         );
addOptional(INPUT,'Anharm'         ,defaultAnharm         );
addOptional(INPUT,'LFreq'          ,defaultLFreq          );
addOptional(INPUT,'L_Index'        ,defaultL_Index        );
addOptional(INPUT,'DiagDisorder'   ,defaultDiagDisorder   );
addOptional(INPUT,'OffDiagDisorder',defaultOffDiagDisorder);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
Vec_1           = INPUT.Results.Vec_1          ;
Vec_2           = INPUT.Results.Vec_2          ;
N_1             = INPUT.Results.N_1            ;
N_2             = INPUT.Results.N_2            ;
NLFreq          = INPUT.Results.NLFreq         ;
Anharm          = INPUT.Results.Anharm         ;
LFreq           = INPUT.Results.LFreq          ;
L_Index         = INPUT.Results.L_Index        ;
DiagDisorder    = INPUT.Results.DiagDisorder   ;
OffDiagDisorder = INPUT.Results.OffDiagDisorder;

%% Extend the monomer along the two directions
S_V1N = SD_TransN(S_Monomer,Vec_1,N_1);
S_V1N_V2M = SD_TransN(S_V1N,Vec_2,N_2);

%% Assign Labeled frequency
S_V1N_V2M.LocFreq         = ones(S_V1N_V2M.Nmodes,1).*NLFreq;
S_V1N_V2M.LocAnharm       = ones(S_V1N_V2M.Nmodes,1).*Anharm;
S_V1N_V2M.DiagDisorder    = ones(S_V1N_V2M.Nmodes,1).*DiagDisorder;
S_V1N_V2M.OffDiagDisorder = ones(S_V1N_V2M.Nmodes,1).*OffDiagDisorder;

if ~ischar(L_Index)
    S_V1N_V2M.LocFreq(L_Index) = LFreq;
end

%% Deal with files name 
S_V1N_V2M.FilesName = [ '2D_Grid_V1' num2str(N_1) '_V2' num2str(N_2)];

%% Save the Translational copy info
Extra.N_Vec1 = N_1;
Extra.N_Vec2 = N_2;
Extra.Vec_1  = Vec_1;
Extra.Vec_2  = Vec_2;

S_V1N_V2M.Extra    = Extra;
S_V1N_V2M.Children = S_Monomer;
