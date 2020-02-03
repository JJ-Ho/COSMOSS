function TCO = ConstructTCO(GUI_Inputs)
% This script construct a Two Coupled Oscillator model using two cabocilic 
% acid molecules. The molecular information was calculated by G09.
% Copyright Jia-Jung Ho, 2014-2020

%% Debug input part
% clear all
% handles = 'Debug';

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultPhi_D1          = 0;
defaultPsi_D1          = 0;
defaultTheta_D1        = 0;
defaultPhi_D2          = 0;
defaultPsi_D2          = 0;
defaultTheta_D2        = 0;
defaultDisplacement    = [0,0,5];
defaultPhi             = 0;
defaultPsi             = 0;
defaultTheta           = 0;

% add Optional inputs / Parameters
addOptional(INPUT,'Phi_D1'         ,defaultPhi_D1         );
addOptional(INPUT,'Psi_D1'         ,defaultPsi_D1         );
addOptional(INPUT,'Theta_D1'       ,defaultTheta_D1       );
addOptional(INPUT,'Phi_D2'         ,defaultPhi_D2         );
addOptional(INPUT,'Psi_D2'         ,defaultPsi_D2         );
addOptional(INPUT,'Theta_D2'       ,defaultTheta_D2       );
addOptional(INPUT,'Displacement'   ,defaultDisplacement   );
addOptional(INPUT,'Phi'            ,defaultPhi            );
addOptional(INPUT,'Psi'            ,defaultPsi            );
addOptional(INPUT,'Theta'          ,defaultTheta          );

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
Phi_D1          = INPUT.Results.Phi_D1;
Psi_D1          = INPUT.Results.Psi_D1;
Theta_D1        = INPUT.Results.Theta_D1;
Phi_D2          = INPUT.Results.Phi_D2;
Psi_D2          = INPUT.Results.Psi_D2;
Theta_D2        = INPUT.Results.Theta_D2;
Displacement    = INPUT.Results.Displacement;
Phi             = INPUT.Results.Phi;
Psi             = INPUT.Results.Psi;
Theta           = INPUT.Results.Theta;

%% Monomer Settings
Num_Modes = 2;

% Monomer XYZ, C, =O, -O, D
XYZ_1 = [0.000,   0.000,   0.000;
         0.000,   0.000,   1.204;
         0.000,   1.142,  -0.730;
         0.000,   1.6380, -0.2197];
         %0.000,   1.905,   0.055];

AtomName = {'C','O','O','H'};             
             
%% Rotate the first chromophore
% Orientation = Orientation/180*pi; % turn to radius unit
Phi_R1   = Phi_D1/180*pi;
Psi_R1   = Psi_D1/180*pi;
Theta_R1 = Theta_D1/180*pi;

Rot_Mat1 = R1_ZYZ_0(Phi_R1,Psi_R1,Theta_R1);
XYZ_1_Rot = (Rot_Mat1*XYZ_1')';

%% Orient the second Acid 
% XYZ_coor = CTCO_1(handles,XYZ_1);

XYZ_2 = XYZ_1;

% Orientation = Orientation/180*pi; % turn to radius unit
Phi_R2   = Phi_D2/180*pi;
Psi_R2   = Psi_D2/180*pi;
Theta_R2 = Theta_D2/180*pi;

Rot_Mat2 = R1_ZYZ_0(Phi_R2,Psi_R2,Theta_R2);
XYZ_2_Rot = (Rot_Mat2*XYZ_2')';
XYZ_2_Rot_Trans = bsxfun(@plus,XYZ_2_Rot,Displacement);

% combine into dimer
SZ = size(XYZ_1);
XYZ_coor = zeros(2,SZ(1),SZ(2));
XYZ_coor(1,:,:) = XYZ_1_Rot;
XYZ_coor(2,:,:) = XYZ_2_Rot_Trans;

XYZ = [ squeeze(XYZ_coor(1,:,:)) ;
        squeeze(XYZ_coor(2,:,:)) ];

AtomName = [AtomName,AtomName]';    
    
%% Define Acid modes coordinate system
Vec_COD      = XYZ_coor(:,2,:)-XYZ_coor(:,1,:);
Vec_COD      = squeeze(Vec_COD);
Vec_COD_Norm = sqrt(sum(abs(Vec_COD).^2,2));
Vec_COD      = bsxfun(@rdivide,Vec_COD,Vec_COD_Norm); % normaliz COD vectors

Vec_COS      = XYZ_coor(:,3,:)-XYZ_coor(:,1,:);
Vec_COS      = squeeze(Vec_COS);
Vec_COS_Norm = sqrt(sum(abs(Vec_COS).^2,2));
Vec_COS      = bsxfun(@rdivide,Vec_COS,Vec_COS_Norm); % normaliz COS vectors

% AmideICenter = squeeze(XYZ_Atom_Rot(:,1,:)) + Vec_COD.*0.665 + Vec_COS.*0.256; % center of amide I mode, ref from Jenny's mathematica code
AcidCenter = squeeze(XYZ_coor(:,1,:));

% Define Lab frame coordinate of each mode
Z_Sim = Vec_COD;
X_Sim = cross(Vec_COS,Z_Sim,2);
X_Sim = bsxfun(@rdivide,X_Sim,sqrt(sum(abs(X_Sim).^2,2))); % normalize
Y_Sim = cross(Z_Sim,X_Sim,2);

XYZ_Sim = [X_Sim(:); Y_Sim(:); Z_Sim(:)]';
XYZ_Sim = reshape(XYZ_Sim,Num_Modes,3,[]); 

%% Calculate the transition dipoles (mu) and Raman tensors (alpha) in the Lab frame
% mu and alpha in G09 simulation frame 
% mu_Mol    = [-5.43155E+00,  1.653500E+01, -3.17990E-05];
% alpha_Mol = [0.395189E-01,  0.239907E+00,  0.000000E+00;
%              0.239907E+00, -0.607441E+00,  0.138181E-05;
%              0.000000E+00,  0.138181E-05,  0.731944E-01 ];

% mu and alpha in [1,0,0] [0,1,0] [0,0,1] frame
% mu_Mol    = [0.0000, -1.8499,  17.3057];
mu_Mol    = [0.0000,  0.0000,  17.4043];
alpha_Mol = [0.0732,  0.0000,   0.0000;
             0.0000,  0.1100,   0.0838;
             0.0000,  0.0838,  -0.6779];
         
mu_Sim    = zeros(Num_Modes,3);
alpha_Sim = zeros(Num_Modes,3,3);

for ii=1:Num_Modes
    
    R_Mol2Sim = squeeze(XYZ_Sim(ii,:,:));
    mu_Sim(ii,:,:)    = R_Mol2Sim * mu_Mol';
    alpha_Sim(ii,:,:) = R_Mol2Sim * alpha_Mol' * R_Mol2Sim';
end    

%% Vectorize Alpha
% alpha_Sim = [N x 3 x3 ]
% for a signle mode, the alpha_Sim: 
% [ XX, XY, XZ ]
% [ YX, YY, YZ ]
% [ ZX, ZY, ZZ ]
% Following the kron convention that the following E J L R beta used, the
% alpha vector need to be in this index order
% [XX,XY, XZ, YX, YY, YZ, ZX, ZY, ZZ]'
% eventhough the Raman tensor we encounter in IR resonance SFG is always 
% symmetric, I am being exta caucious here to make the indexing right.
alpha = reshape(permute(alpha_Sim,[1,3,2]),[Num_Modes,9]);

%% Constructe the StructureData Object with given Structural information
TCO = StructureData;
TCO.XYZ             = XYZ;
TCO.AtomName        = AtomName;
TCO.LocCenter       = AcidCenter;
TCO.LocMu           = mu_Sim;
TCO.LocAlpha        = alpha; % raman tensor vector form [N x 9]

%% Post process of the dimer
% Rotate the molecule
Phi_r   = Phi/180*pi;
Psi_r   = Psi/180*pi;
Theta_r = Theta/180*pi;
R1      = R1_ZYZ_0(Phi_r,Psi_r,Theta_r);
TCO     = SD_Rot(TCO,R1);
            
% Move the center of the two Carbon atoms to the  origin (0,0,0)
CV  = mean(TCO.XYZ([1,5],:));
TCO = SD_Trans(TCO,-CV);

%% Add Model dependent properties
TCO.FilesName  = 'Acid Dimer';
TCO.GUI_Inputs = GUI_Inputs;
TCO.hPlotFunc  = @PlotXYZfiles_TCO;
TCO.GUI_Inputs = GUI_Inputs;

%% Calculate One Exciton Hamiltoian
TCO = SD_1ExH(TCO); 

