function NCO = ConstructNCO(GUI_Inputs)
% This script construct a Two Coupled Oscillator model using two cabocilic 
% acid molecules. The molecular information was calculated by G09.
% Copyright Jia-Jung Ho, 2014-2020

%% Debug input part
% clear all
% GUI_Inputs.mode = 'Debug';
% S = Data_TCO.GUI_Inputs.StructureInputs;

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultPhi          = 0;
defaultPsi          = 0;
defaultTheta        = 0;
defaultLocFreqType  = 'NCO';

% add Optional inputs / Parameters
addOptional(INPUT,'Phi'         ,defaultPhi         );
addOptional(INPUT,'Psi'         ,defaultPsi         );
addOptional(INPUT,'Theta'       ,defaultTheta       );
addOptional(INPUT,'LocFreqType' ,defaultLocFreqType );

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
Phi          = INPUT.Results.Phi;
Psi          = INPUT.Results.Psi;
Theta        = INPUT.Results.Theta;
LocFreqType  = INPUT.Results.LocFreqType;

%% Parse the structure info from the GUI input table
S = GUI_Inputs.StructureInputs;
N_Modes_toggle = logical(S(:,1));
S = S(N_Modes_toggle,:);

Phi_monomer   = S(:,2)./180*pi; % turn to radius unit
Psi_monomer   = S(:,3)./180*pi; % turn to radius unit
Theta_monomer = S(:,4)./180*pi; % turn to radius unit
Displacement  = S(:,5:7);
ScalingFactor = S(:,8);

Num_Modes = size(Phi_monomer,1);

%% Monomer Settings, take carboxylic acid group as example
% Monomer XYZ, C, =O, -O, D
XYZ_monomer = [0.000,   0.000,   0.602;
               0.000,   0.000,  -0.602];

AtomName = {'N','N'}';   

AcidCenter = XYZ_monomer(1,:);

mu_Mol    = [0.0000,  0.0000,  17.4043];
alpha_Mol = [0.0732,  0.0000,   0.0000, 0.0000,  0.1100,   0.0838, 0.0000,  0.0838,  -0.6779];

%% Create the N-coupled oscillator system
NCO = StructureData;
NCO.XYZ       = XYZ_monomer;
NCO.AtomName  = AtomName;
NCO.LocCenter = NCO.CoM;
NCO.LocMu     = mu_Mol   .*ScalingFactor(1);
NCO.LocAlpha  = alpha_Mol.*ScalingFactor(1);
NCO.Beta      = 0;

Rot_M = R1_ZYZ_0(Phi_monomer(1,:),Psi_monomer(1,:),Theta_monomer(1,:));
NCO = SD_Rot(NCO,Rot_M);

NCO = SD_Trans(NCO,Displacement(1,:));

for i = 2:Num_Modes
    Monomer_tmp = StructureData;
    Monomer_tmp.XYZ       = XYZ_monomer;
    Monomer_tmp.AtomName  = AtomName;
    Monomer_tmp.LocCenter = AcidCenter;
    Monomer_tmp.LocMu     = mu_Mol   .*ScalingFactor(i);
    Monomer_tmp.LocAlpha  = alpha_Mol.*ScalingFactor(i);
    Monomer_tmp.Beta      = 0;
    
    Rot_M = R1_ZYZ_0(Phi_monomer(i,:),Psi_monomer(i,:),Theta_monomer(i,:));
    Monomer_tmp = SD_Rot(Monomer_tmp,Rot_M);
    
    Monomer_tmp = SD_Trans(Monomer_tmp,Displacement(i,:));
    
    NCO = SD_Comb2(NCO,Monomer_tmp,'Zero',0); % give a dummy coupling model, will deal with it later
end

% remove unnecessary Children
NCO.Children = [];

%% Post process of the dimer
% Rotate the molecule
Phi_r   = Phi/180*pi;
Psi_r   = Psi/180*pi;
Theta_r = Theta/180*pi;
R1      = R1_ZYZ_0(Phi_r,Psi_r,Theta_r);
NCO     = SD_Rot(NCO,R1);
            
%% Feeds the default GUI_Options back into the StructureData
GUI_Inputs.LocFreqType = LocFreqType;

%% Add Model dependent properties
NCO.FilesName  = 'Acid N-mer';
NCO.GUI_Inputs = GUI_Inputs;
NCO.hPlotFunc  = @PlotXYZ_Grid;
NCO.GUI_Inputs = GUI_Inputs;

%% Calculate One Exciton Hamiltonian
NCO = SD_1ExH(NCO); 

