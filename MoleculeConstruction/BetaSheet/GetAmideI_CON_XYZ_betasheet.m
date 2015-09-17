function Output = GetAmideI_CON_XYZ_betasheet(XYZ,varargin)

%% GetAmideI
% Output = GetAmideI(Betasheet_XYZ,[Isotopes,LabelFreq])
%  
% This Script read pdb file and isolate the coordinates of (C,O,N) atoms of
% the amide I mode then output the [center axes] of each local amide I
% mode.
% 
% PDB_Code: pdb code of molecule with out extension (.pdb)
% Isotopes: array of isotoped A.A. position
% LabelFreq: Labeled frequencies, same length of Isotopes
% 

% Todo: integrate InputParser
% 
% 
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.3  140122  Change index of XYZ according to the modification of
%                   BetaSheet generation code.
% 
% Ver. 1.2  140119  Change local mode freq from 1645 to 1644
% 
% Ver. 1.1  140110  Add "NumLocMode" as output
%                   Capitalize the first letter output variables
% 
% Ver. 1.0  140102  Modified from GetAmideI_CON.m to extract model
%                   betasheet
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014


ModeNum = size(XYZ,1);

Vec_CO = XYZ(:,:,2)-XYZ(:,:,1);
Vec_CO = squeeze(Vec_CO);
Vec_CO_Norm = sqrt(sum(abs(Vec_CO).^2,2));
Vec_CO = bsxfun(@rdivide,Vec_CO,Vec_CO_Norm); % normaliz CO vectors

Vec_CN = XYZ(:,:,3)-XYZ(:,:,2);
Vec_CN = squeeze(Vec_CN);
Vec_CN_Norm = sqrt(sum(abs(Vec_CN).^2,2));
Vec_CN = bsxfun(@rdivide,Vec_CN,Vec_CN_Norm); % normaliz CN vectors

AmideICenter = squeeze(XYZ(:,:,1)) + Vec_CO.*0.665 + Vec_CN.*0.256; % center of amide I mode, ref from Jenny's mathematica code
%% Define Lab frame coordinate of each mode

Z_Sim = Vec_CO;
X_Sim = cross(Z_Sim,Vec_CN,2);
X_Sim = bsxfun(@rdivide,X_Sim,sqrt(sum(abs(X_Sim).^2,2))); % normalize
Y_Sim = cross(X_Sim,Z_Sim,2);

XYZ_Sim = [X_Sim(:); Y_Sim(:); Z_Sim(:)]';
XYZ_Sim = reshape(XYZ_Sim,ModeNum,3,[]); 

% Lab frame output format:
% 
% [ mode1_X_vec, mode1_Y_vec, mode1_Z_vec;
%   mode2_X_vec, mode2_Y_vec, mode2_Z_vec;
%                                     ...]
% size(XYZ_Lab) = Mode Number * Number of Atom in the mode * 3

%% Calculate rotational matrixes to rotate Molecular frame into Lab frame for each modes
% Actually I don't even need to calculate! since:
% XYZ_Lab = Rot*xyz_Mol = Rot*Identity = Rot

%% Calculate the transition dipole (mu) and Raman tensor (alpha) in Lab frame

% Transition dipole Angle
% mu_Mol_angle = 27.5*pi/180; from Jenny's paper (10.1021/jp408064b)
mu_Mol_angle = 20*pi/180; % from D. Strasfeld's paper (10.1021/jp9072203)
% mu_Mol_angle = 25*pi/180; % from Cho's paper

mu_Mol = [0, sin(mu_Mol_angle), -cos(mu_Mol_angle)].* 16.1016; % this scale factor was calculated from DFT simulation, check 140207_NMA


alpha_angle = 34*pi/180;
alpha = [ 0.05,0,0; 0,0.2,0; 0,0,1].*5;
Rot_to_Mol = [ 1,0,0; 0,cos(alpha_angle),-sin(alpha_angle); 0, sin(alpha_angle),cos(alpha_angle) ];
alpha_Mol = Rot_to_Mol*alpha*Rot_to_Mol';

mu_Sim = zeros(ModeNum,3);
alpha_Sim = zeros(ModeNum,3,3);

for i=1:ModeNum
    mu_Sim(i,:,:) = squeeze(XYZ_Sim(i,:,:))*mu_Mol';
    alpha_Sim(i,:,:) = squeeze(XYZ_Sim(i,:,:))*alpha_Mol'*squeeze(XYZ_Sim(i,:,:))'; % see tensor rotation in evernote
end

% Reduce Raman tensor from 9 elements to 6 unique elements (since it is symmetic)
% --------------------------------------------------------------------------------------- 
%                       [ XX XY XZ ]        
% alpha_Sim(i,:,:) =    [ YX YY YZ ]    => alpha_Sim_reduced(i,:) = [XX YY ZZ XY YZ XZ]
%                       [ ZX ZY ZZ ]
% ---------------------------------------------------------------------------------------

alpha_Sim_reduced = zeros(ModeNum,6);
for j=1:ModeNum
    alpha_Sim_reduced(j,1:3) = diag(squeeze(alpha_Sim(j,:,:)),0);
    alpha_Sim_reduced(j,4:5) = diag(squeeze(alpha_Sim(j,:,:)),1);
    alpha_Sim_reduced(j,6)   = diag(squeeze(alpha_Sim(j,:,:)),2);
end
    

%% Mode frequency and anharmonicity
% varargin{1},varargin{2} = Isotopes,LabelFreq

AmideIFreq = ones(ModeNum,1)*1644; % 1644 for unlabeled beta sheet
% AmideIFreq = ones(ModeNum,1)*1665.8; % from Cho's paper

if nargin == 3
    AmideIFreq(varargin{1}) = varargin{2};
end

% anharmonicity = 14, for betasheet, ref from Lauren's code
AmideIAnharm = ones(ModeNum,1)*14;


%% Output Structure
Output.XYZ              = XYZ;
Output.center           = AmideICenter;
Output.freq             = AmideIFreq;
Output.Num_Modes        = ModeNum;
Output.anharm           = AmideIAnharm;
Output.mu               = mu_Sim;
Output.alpha_matrix     = alpha_Sim;
Output.alpha            = reshape(alpha_Sim,[ModeNum,9]); % non-reduced alpha "vector"
% Output.Alpha_Reduced    = alpha_Sim_reduced;
