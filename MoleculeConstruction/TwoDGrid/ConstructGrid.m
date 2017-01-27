function Output = ConstructGrid(Monomer,GUI_Inputs)
% ConstructGrid
%
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.1  150426  Copy from 2DSFG_Ester project
%                   Update debug part
%                   Implement inputparser
% 
% Ver. 1.0  140903  Modify from GetAcid.m 
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014-2016

%% Debug input part
% clear all
% GUI_Inputs.debug = 1;

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultDelta_Phi   = 0;
defaultDelta_Psi   = 0;
defaultDelta_Theta = 0;
defaultVec_1       = [7,0,0];
defaultVec_2       = [0,4,0];
defaultN_1         = 2;
defaultN_2         = 3;
defaultNLFreq      = 1720;
defaultAnharm      = 20;
defaultLFreq       = 1700;
defaultL_Index     = 'None';
defaultMute_Ind    = 'Non';

% add Optional inputs / Parameters
addOptional(INPUT,'Delta_Phi'  ,defaultDelta_Phi  );
addOptional(INPUT,'Delta_Psi'  ,defaultDelta_Psi  );
addOptional(INPUT,'Delta_Theta',defaultDelta_Theta);
addOptional(INPUT,'Vec_1'      ,defaultVec_1      );
addOptional(INPUT,'Vec_2'      ,defaultVec_2      );
addOptional(INPUT,'N_1'        ,defaultN_1        );
addOptional(INPUT,'N_2'        ,defaultN_2        );
addOptional(INPUT,'NLFreq'     ,defaultNLFreq     );
addOptional(INPUT,'Anharm'     ,defaultAnharm     );
addOptional(INPUT,'LFreq'      ,defaultLFreq      );
addOptional(INPUT,'L_Index'    ,defaultL_Index    );
addOptional(INPUT,'Mute_Ind'   ,defaultMute_Ind   );

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
Delta_Phi   = INPUT.Results.Delta_Phi  ;
Delta_Psi   = INPUT.Results.Delta_Psi  ;
Delta_Theta = INPUT.Results.Delta_Theta;
Vec_1       = INPUT.Results.Vec_1      ;
Vec_2       = INPUT.Results.Vec_2      ;
N_1         = INPUT.Results.N_1        ;
N_2         = INPUT.Results.N_2        ;
NLFreq      = INPUT.Results.NLFreq     ;
Anharm      = INPUT.Results.Anharm     ;
LFreq       = INPUT.Results.LFreq      ;
L_Index     = INPUT.Results.L_Index    ;
Mute_Ind    = INPUT.Results.Mute_Ind   ;

%% Read G09 structure and reassign variables

XYZ_G09       = Monomer.XYZ;
Atom_Num_G09  = Monomer.NAtoms;
mu_Mol_G09    = Monomer.LocMu; % size [N x 3]
alpha_Mol_G09 = Monomer.LocAlpha; % note: RamanV = [N x 9], index: [xx xy xz yx yy yz zx zy zz]
Atom_Name_G09 = Monomer.AtomName;
Center_Ind    = Monomer.Extra.Center_Ind;

%% Define usefule constants

% define 2D grid vector
% Vec_1 = a.*[1,0,0];
% Vec_2 = b.*[0,1,0];

Num_Modes = N_1*N_2;

Phi_Fluc   =   Delta_Phi*randn(Num_Modes,1)./180*pi;
Psi_Fluc   =   Delta_Psi*randn(Num_Modes,1)./180*pi;
Theta_Fluc = Delta_Theta*randn(Num_Modes,1)./180*pi;

%% Creating Translation Copy
XYZ_Grid_M = zeros(Atom_Num_G09,3,N_1,N_2);

mu_Sim    = zeros(Num_Modes,3);
alpha_Sim = zeros(Num_Modes,9);

for j = 1:N_2
    for i = 1:N_1

        if any(Mute_Ind == i+(j-1)*N_1)
            % Mute some molecule in the grid
            XYZ_Grid_M(:,:,i,j)        = zeros(Atom_Num_G09,3);
            mu_Sim(i+(j-1)*N_1,:)    = zeros(3,1);
            alpha_Sim(i+(j-1)*N_1,:) = zeros(9,1);
        else     
            % Add random orientation fluctuation
            Phi_Ind = Phi_Fluc(i+(j-1)*N_1);
            Psi_Ind = Psi_Fluc(i+(j-1)*N_1);
            Theta_Ind = Theta_Fluc(i+(j-1)*N_1);
            R1_Fluc = R1_ZYZ_0(Phi_Ind,Psi_Ind,Theta_Ind);
            R2_Fluc = R2_ZYZ_0(Phi_Ind,Psi_Ind,Theta_Ind);

            % XYZ
            XYZ_R  = (R1_Fluc*XYZ_G09')'; % apply random fluctuation 
            TransV = (i-1)*Vec_1 + (j-1)*Vec_2;
            XYZ_Grid_M(:,:,i,j) = bsxfun(@plus,XYZ_R,TransV);

            % Mu & Alpha
               mu_Sim(i+(j-1)*N_1,:) = (R1_Fluc * mu_Mol_G09'    )'; % note    mu_Mol_G09 = [N_Mode*3]
            alpha_Sim(i+(j-1)*N_1,:) = (R2_Fluc * alpha_Mol_G09' )'; % note alpha_Mol_G09 = [N_Mode*9]
        end
    end
end

XYZ_Grid = reshape(permute(XYZ_Grid_M,[1,3,4,2]),[],3);

%% Extend Atom_Name
Atom_Name = repmat(Atom_Name_G09,Num_Modes,1);

%% Create Translational copy of Center

Center_M = sum(XYZ_Grid_M(Center_Ind,:,:,:),1);
Center = reshape(permute(Center_M,[1,3,4,2]),[],3);

%% Define Mode frequency and anharmonicity

Loc_Freq = ones(Num_Modes,1)* NLFreq;

if ~ischar(L_Index)
    Loc_Freq(L_Index) = LFreq.*ones(size(L_Index));
end

% Anharm = ones(Num_Modes,1)*(2*Freq_G09.Fundamental - Freq_G09.Overtone);
Anharm = Anharm.*ones(Num_Modes,1);

%% AtomSerNo

% AtomSerNo = zeros(Num_Modes,3);
% for jj = 1:Num_Modes
%     Shift_Ind = (jj-1)*Atom_Num_G09;
%     AtomSerNo(jj,:) = [7+Shift_Ind,10+Shift_Ind,8+Shift_Ind];
% end

%% Deal with files name 
Grid_FilesName = [ '2D_Grid_V1' num2str(N_1) 'V2' num2str(N_2)];

%% Output Structure
Output = StructureData;

Output.XYZ       = XYZ_Grid;
Output.AtomName  = Atom_Name;
Output.COM       = sum(XYZ_Grid,1)./size(XYZ_Grid,1);

Output.LocCenter = Center;
Output.LocFreq   = Loc_Freq;
Output.LocAnharm = Anharm;
Output.LocMu     = mu_Sim; % size [N x 3]
Output.LocAlpha  = alpha_Sim; % raman tensor vector form [N x 9]

Output.FilesName = Grid_FilesName;

Extra.N_Vec1 = N_1;
Extra.N_Vec2 = N_2;
Extra.Vec_1  = Vec_1;
Extra.Vec_2  = Vec_2;
Output.Extra = Extra;

% Output.center       = Center;
% Output.freq         = Loc_Freq;
% Output.anharm       = Anharm;
% Output.mu           = mu_Sim; % size [N x 3]
% Output.alpha        = alpha_Sim; % note: RamanV = [N x 9], index: [xx xy xz yx yy yz zx zy zz]
% Output.alpha_matrix = reshape(alpha_Sim,[Num_Modes,3,3]);
% Output.AtomSerNo    = AtomSerNo;
% Output.Num_Modes    = Num_Modes;
% Output.XYZ          = XYZ_Grid;
% Output.FilesName    = Grid_FilesName;
% Output.Monomer      = Monomer;
% Output.N_Vec1       = N_1;
% Output.N_Vec2       = N_2;
% Output.Vec_1        = Vec_1;
% Output.Vec_2        = Vec_2;
% Output.Atom_Name    = Atom_Name;
