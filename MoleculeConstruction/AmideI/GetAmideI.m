function Output = GetAmideI(Num_Atoms,XYZ,AtomName,FilesName,varargin)

%% GetAmideI
% Output = GetAmideI([Isotopes,LabelFreq])
%  
% This Script read pdb file and isolate the coordinates of (C,O,N) atoms of
% the amide I mode then output the [center axes] of each local amide I
% mode.
% 
% Isotopes: array of isotoped A.A. position
% LabelFreq: Labeled frequencies, same length of Isotopes
% 
% If ther's no isotope labeling, just need single PDB_code input
% 
% Example:
% 
%   Output = GetAmideI(2RRI);
% 
%   Output = GetAmideI(2RRI,[10,13],[1650,1643]);
% 

% ------- Version log -----------------------------------------------------
% 
% Ver. 3.1  141021  Copy from SFG_AmideI_GUI
%                   Add Input parser
% 
% Ver. 3.0  140825  Correct the TDV direction from pointing to N stom to
%                   pointing away from N atom. See DFT calculation of 
%                   140207_NMA for TDV vector direction. 
% 
% Ver. 2.9  140824  Add "Dimenesion" to determine sum/cross dimention for 
%                   one mode case.
% 
% Ver. 2.8  140822  Clean up the code accessing tedious cell structure
%                   Correct Aminde I modes coordinate system and the TDV
%                   direction.
% 
% Ver. 2.7  140605  Add XYZ as output
% 
% Ver. 2.6  140603  Add uigetfile to select input; Adjust output names
% 
% Ver. 2.5  131108  change name of ouput variable "mu_orig" to "mu"
% 
% Ver. 2.4  130925  fix home filder dependancy in "pdb file location" part
%                   export alpha in both matrix form and vector form.
% 
% Ver. 2.3  130813  AmideI miss assignmant fixed by using while loop
% 
% Ver. 2.2  130731  use varargin to avoid error when there's no labeling
%                   inputs. varargin{1} = label index array; varargin{2} = 
%                   label frequency array.
% 
% Ver. 2.1  130723  Reduce Raman tensor to 6 unique element after final
%                   rotation. The output of unique Raman tensor is
%                   linearized 6*1 vector.
%                   And rename all "Lab" frames variables to "Sim" frame
%                   since the coordinate transform is from Molecule frame
%                   to Simulation frame here. 
% 
% Ver. 2.0  130706  Change this script to function, move plot molecular
%                   part to upper level => TwoDSFG_Simulation. And remove
%                   plot_toggle.
% 
% Ver. 1.1  130619  Atom_Data = TT.Model.Atom
% 
% Ver. 1.0  130604 
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

%% Debug input part
% clear all
% handles = 'Debug';

% if ~isstruct(handles)
%     Phi_D = 0;
%     Psi_D = 0;
%     Theta_D = 0;
% end

%% Inputs parser
INPUT = inputParser;

% Default values
defaultPhi_D        = 0;
defaultPsi_D        = 0;
defaultTheta_D      = 0;
defaultNLFreq       = 1644;
defaultAnharm       = 12;

% add Optional inputs / Parameters
addOptional(INPUT,'Phi_D'       ,defaultPhi_D       );
addOptional(INPUT,'Psi_D'       ,defaultPsi_D       );
addOptional(INPUT,'Theta_D'     ,defaultTheta_D     );
addOptional(INPUT,'NLFreq',defaultNLFreq);
addOptional(INPUT,'Anharm',defaultAnharm);

parse(INPUT,varargin{:});

% Re-assign variable names
Phi_D         = INPUT.Results.Phi_D;
Psi_D         = INPUT.Results.Psi_D;
Theta_D       = INPUT.Results.Theta_D;
NLFreq        = INPUT.Results.NLFreq;
Anharm        = INPUT.Results.Anharm;


%% Get pdb file location 
% PWD = pwd;
% PDB_Path = [PWD, '/PDB_files/'];
% 
% [FilesName,PathName,FilterIndex] = uigetfile({'*.pdb','PDB file'; ...
%                                               '*,*','All Files'},...
%                                              'Select inputs',PDB_Path);

%% Parse molecule structure

% switch FilterIndex
%     case 1 % for PDB files
%         PDB = pdbread([PathName FilesName]);
%         Atom_Data = PDB.Model.Atom;
%         Num_Atoms = size(Atom_Data,2);
% 
%         % Get coordination data
%         XYZ = zeros(Num_Atoms,3);
%         AtomName = cell(Num_Atoms,1);
%         for II = 1:Num_Atoms
%             A = Atom_Data(II);
%             XYZ(II,:) = [A.X, A.Y, A.Z];
%             AtomName{II} = Atom_Data(II).AtomName;
%         end

% I decided not to read XYZ file diretly since it is hard to identify amide
% mode without elemnt name information, such as "CA, HD1,..."  
%     case 2 % for XYZ files
%         fid=fopen([PathName FilesName]);
%         Raw = textscan(fid,'%s%f64%f64%f64','HeaderLines',2);
%         fclose(fid);
%         
%         XYZ =[Raw{2},Raw{3},Raw{4}];
%         AtomName = Raw{1};
%         Num_Atoms = size(XYZ,1);
% end
        


%% find corresponding atoms for Amide I mode

AmideI_C_AtomSerNo = [];

for ii=1:Num_Atoms
    if strcmp(AtomName{ii},'C')
        AmideI_C_AtomSerNo = [AmideI_C_AtomSerNo;ii];
    end
end

Num_AmideC = length(AmideI_C_AtomSerNo);
AmideI_O_AtomSerNo = zeros(Num_AmideC,1);
AmideI_N_AtomSerNo = zeros(Num_AmideC,1);

Count = 1;
for jj=1:length(AmideI_C_AtomSerNo)
    
    FoundAmideO = 'n';
    AtomIndex = AmideI_C_AtomSerNo(jj)+1;
    
    while strcmp(FoundAmideO,'n');
        if strcmp(AtomName{AtomIndex},'O')
            AmideI_O_AtomSerNo(Count) = AtomIndex;
            
            FoundAmideO = 'y';
            Count = Count + 1;
        else
            AtomIndex = AtomIndex + 1;
            if AtomIndex > Num_Atoms
                break
            end
        end 
    end
end

Count = 1;
for kk=1:length(AmideI_C_AtomSerNo)
    
    FoundAmideN = 'n';
    AtomIndex = AmideI_O_AtomSerNo(kk)+1;
    
    while strcmp(FoundAmideN,'n');
        if strcmp(AtomName{AtomIndex},'N')
            AmideI_N_AtomSerNo(Count) = AtomIndex;
            
            FoundAmideN = 'y';
            Count = Count + 1;
        else
            AtomIndex = AtomIndex + 1;
            if AtomIndex > Num_Atoms
                break
            end
        end 
    end
end
            
        
AmideIAtomSerNo = [AmideI_C_AtomSerNo';AmideI_O_AtomSerNo';AmideI_N_AtomSerNo']';
% delete rows without N atom presenting, this is always the last line
AmideIAtomSerNo = AmideIAtomSerNo(logical(AmideIAtomSerNo(:,3)),:);

%% Read molecule XYZ and rotate molecule in Molecule frame
Num_Modes = size(AmideIAtomSerNo,1);

XYZ_Orig = XYZ(AmideIAtomSerNo(:),:);

% Orientation = Orientation/180*pi; % turn to radius unit
Phi_R   = Phi_D/180*pi;
Psi_R   = Psi_D/180*pi;
Theta_R = Theta_D/180*pi;

Rot_Mat = R1_ZYZ_0(Phi_R,Psi_R,Theta_R);
XYZ_Atom_Rot = (Rot_Mat*XYZ_Orig')';


XYZ_Atom_Rot = reshape(XYZ_Atom_Rot,Num_Modes,3,[]); % [C_xyz O_xyz N_xyz ; C_xyz, O_xyz, N_xyz ; ...]

%% Define Aminde I modes coordinate system
switch Num_Modes
    case 1
        Dimenesion = 1;
    otherwise
        Dimenesion = 2;
end

Vec_CO = XYZ_Atom_Rot(:,2,:)-XYZ_Atom_Rot(:,1,:);
Vec_CO = squeeze(Vec_CO);
Vec_CO_Norm = sqrt(sum(abs(Vec_CO).^2,Dimenesion));
Vec_CO = bsxfun(@rdivide,Vec_CO,Vec_CO_Norm); % normaliz CO vectors

Vec_CN = XYZ_Atom_Rot(:,3,:)-XYZ_Atom_Rot(:,1,:);
Vec_CN= squeeze(Vec_CN);
Vec_CN_Norm = sqrt(sum(abs(Vec_CN).^2,Dimenesion));
Vec_CN = bsxfun(@rdivide,Vec_CN,Vec_CN_Norm); % normaliz CN vectors

AmideICenter = squeeze(XYZ_Atom_Rot(:,1,:)) + Vec_CO.*0.665 + Vec_CN.*0.256; % center of amide I mode, ref from Jenny's mathematica code

% Define Lab frame coordinate of each mode
Z_Sim = Vec_CO;
X_Sim = cross(Z_Sim,Vec_CN,Dimenesion);
X_Sim = bsxfun(@rdivide,X_Sim,sqrt(sum(abs(X_Sim).^2,Dimenesion))); % normalize
Y_Sim = cross(Z_Sim,X_Sim,Dimenesion);

XYZ_Sim = [X_Sim(:); Y_Sim(:); Z_Sim(:)]';
XYZ_Sim = reshape(XYZ_Sim,Num_Modes,3,[]); 


%% Calculate the transition dipoles (mu) and Raman tensors (alpha) in the Lab frame
mu_Mol_angle = 27.5*pi/180;
mu_Mol = [0, sin(mu_Mol_angle), cos(mu_Mol_angle)].* 16.1016; % this scale factor was calculated from DFT simulation, check 140207_NMA


alpha_angle = 34*pi/180;
alpha = [ 0.05,0,0; 0,0.2,0; 0,0,1].*5;
Rot_to_Mol = [ 1,0,0; 0,cos(alpha_angle),-sin(alpha_angle); 0, sin(alpha_angle),cos(alpha_angle) ];
alpha_Mol = Rot_to_Mol*alpha*Rot_to_Mol';

mu_Sim    = zeros(Num_Modes,3);
alpha_Sim = zeros(Num_Modes,3,3);

for ii=1:Num_Modes
    mu_Sim(ii,:,:)    = squeeze(XYZ_Sim(ii,:,:))*mu_Mol';
    alpha_Sim(ii,:,:) = squeeze(XYZ_Sim(ii,:,:))*alpha_Mol'*squeeze(XYZ_Sim(ii,:,:))';
end    

%% Define Mode frequency and anharmonicity

AmideIFreq = ones(Num_Modes,1).*NLFreq;

% anharmonicity = 12
AmideIAnharm = ones(Num_Modes,1).*Anharm;

%% Output Structure
Output.center       = AmideICenter;
Output.freq         = AmideIFreq;
Output.anharm       = AmideIAnharm;
Output.mu           = mu_Sim;
Output.alpha        = reshape(alpha_Sim,[Num_Modes,9]); % non-reduced alpha "vector"
Output.alpha_matrix = alpha_Sim;
Output.AtomSerNo    = AmideIAtomSerNo;
Output.Num_Modes    = Num_Modes;
Output.XYZ          = XYZ;
Output.FilesName    = FilesName;


