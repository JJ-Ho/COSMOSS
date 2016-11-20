function Output = GetAmideI(Num_Atoms,XYZ,AtomName,FilesName,GUI_Inputs)
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
% varargin = {'Phi_D',0};
% 
% SheetType = 'Anti';
% N_Residue= 7;
% N_Strand = 5;
% TransV = [0,0,4.75];
% TwistV = [0,0,4];
% BB = ConstuctBetaSheet(SheetType,N_Residue,N_Strand,TransV,TwistV);
% 
% Num_Atoms = BB.Num_Atoms;
% XYZ       = BB.XYZ;
% AtomName  = BB.AtomName;
% FilesName = BB.FilesName;

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultPhi_D        = 0;
defaultPsi_D        = 0;
defaultTheta_D      = 0;
defaultNLFreq       = 1644;
defaultAnharm       = 12;
defaultLFreq        = 1604;
defaultL_Index      = 'None';

% add Optional inputs / Parameters
addOptional(INPUT,'Phi_D'       ,defaultPhi_D       );
addOptional(INPUT,'Psi_D'       ,defaultPsi_D       );
addOptional(INPUT,'Theta_D'     ,defaultTheta_D     );
addOptional(INPUT,'NLFreq'      ,defaultNLFreq      );
addOptional(INPUT,'Anharm'      ,defaultAnharm      );
addOptional(INPUT,'LFreq'       ,defaultLFreq       );
addOptional(INPUT,'L_Index'     ,defaultL_Index     );

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
Phi_D         = INPUT.Results.Phi_D;
Psi_D         = INPUT.Results.Psi_D;
Theta_D       = INPUT.Results.Theta_D;
NLFreq        = INPUT.Results.NLFreq;
Anharm        = INPUT.Results.Anharm;
LFreq         = INPUT.Results.LFreq;
L_Index       = INPUT.Results.L_Index;

%% find corresponding atoms for Amide I mode
% The peptide atom patter in pdb file is -[N-CA-C-O-(side chain atoms)]-
% This mode selection algorithm start by searching C atoms and fllowing
% with the search of O,N,CA before hitting next C atom.

% find C atoms
AmideI_C_AtomSerNo = [];
for ii=1:Num_Atoms
    if strcmp(AtomName{ii},'C')
        AmideI_C_AtomSerNo = [AmideI_C_AtomSerNo;ii];
    end
end


Num_AmideC = length(AmideI_C_AtomSerNo);
AmideI_O_AtomSerNo  = nan(Num_AmideC,1);
AmideI_N_AtomSerNo  = nan(Num_AmideC,1);
AmideI_CA_AtomSerNo = nan(Num_AmideC,1);

% find O atom within two the C atoms segment
Count = 1;
for jj=1:Num_AmideC
    
    FoundAmideO = 'n';
    AtomIndex = AmideI_C_AtomSerNo(jj)+1;
    
    while strcmp(FoundAmideO,'n');
        
        % determine if the segament between C(jj) and C(jj+1) is finished   
        switch jj
            case Num_AmideC             
                if AtomIndex > Num_Atoms
                    Count = Count + 1;
                    break
                end
            otherwise
                if AtomIndex > AmideI_C_AtomSerNo(jj+1)
                    Count = Count + 1;
                    break
                end
        end    
        % if the segment is not finished yet continue the search
        if strcmp(AtomName{AtomIndex},'O')
            AmideI_O_AtomSerNo(Count) = AtomIndex;
            
            FoundAmideO = 'y';
            Count = Count + 1;
        else
            AtomIndex = AtomIndex + 1;
        end
    end
end

% find N atom within two the C atoms segment
Count = 1;
for kk=1:Num_AmideC
    
    FoundAmideN = 'n';
    AtomIndex = AmideI_O_AtomSerNo(kk)+1;
    
    while strcmp(FoundAmideN,'n');
        % determine if the segament between C(jj) and C(jj+1) is finished
        switch kk
            case Num_AmideC              
                if AtomIndex > Num_Atoms
                    Count = Count + 1;
                    break
                end
            otherwise
                if AtomIndex > AmideI_C_AtomSerNo(kk+1)
                    Count = Count + 1;
                    break
                end
        end
        % if the segment is not finished yet continue the search
        if strcmp(AtomName{AtomIndex},'N')
            AmideI_N_AtomSerNo(Count) = AtomIndex;
            
            FoundAmideN = 'y';
            Count = Count + 1;
        else
            AtomIndex = AtomIndex + 1;
        end 
    end
end


% find CA atom within two the C atoms segment
Count = 1;
for ll=1:Num_AmideC
    
    FoundAmideCA = 'n';
    AtomIndex = AmideI_N_AtomSerNo(ll)+1;
    
    while strcmp(FoundAmideCA,'n');
        % break while loop if the previous Atom is not exixt
        if isnan(AtomIndex)
            Count = Count + 1;
            break
        end
        
        % determine if the segament between C(jj) and C(jj+1) is finished
        switch ll
            case Num_AmideC              
                if AtomIndex > Num_Atoms
                    Count = Count + 1;
                    break
                end
            otherwise
                if AtomIndex > AmideI_C_AtomSerNo(ll+1)
                    AmideI_CA_AtomSerNo(Count) = nan;
                    Count = Count + 1;
                    break
                end
        end
        % if the segment is not finished yet continue the search
        if strcmp(AtomName{AtomIndex},'CA')
            AmideI_CA_AtomSerNo(Count) = AtomIndex;
            
            FoundAmideCA = 'y';
            Count = Count + 1;
        else
            AtomIndex = AtomIndex + 1;
        end 
    end
end
       
AmideIAtomSerNo = [AmideI_C_AtomSerNo';AmideI_O_AtomSerNo';AmideI_N_AtomSerNo';AmideI_CA_AtomSerNo']';

% delete NaN rows to eliminate non-amide modes
AmideIAtomSerNo(any(isnan(AmideIAtomSerNo),2),:) = [];

%% Read molecule XYZ and rotate molecule in Molecule frame
Num_Modes = size(AmideIAtomSerNo,1);

XYZ_Orig = XYZ(AmideIAtomSerNo(:),:);
Mol_Frame_Orig = [1,0,0;0,1,0;0,0,1];

% Orientation = Orientation/180*pi; % turn to radius unit
Phi_R   = Phi_D/180*pi;
Psi_R   = Psi_D/180*pi;
Theta_R = Theta_D/180*pi;

Rot_Mat = R1_ZYZ_0(Phi_R,Psi_R,Theta_R);
XYZ_Atom_Rot = (Rot_Mat*XYZ_Orig')';

XYZ_Atom_Rot = reshape(XYZ_Atom_Rot,Num_Modes,[],3); % [C_xyz O_xyz N_xyz CA_xyz; C_xyz, O_xyz, N_xyz, CA_xyz ; ...]

% rotate the whole XYZ coordinate
XYZ_Rot = (Rot_Mat*XYZ')';
Mol_Frame_Rot = (Rot_Mat\Mol_Frame_Orig')'; % axis rotation is inv(Rot_M)*Axis_M = Rot_M\Axis_M

%% Define Aminde I modes coordinate system
switch Num_Modes
    case 1
        Dimenesion = 1;
    otherwise
        Dimenesion = 2;
end

Vec_CO = XYZ_Atom_Rot(:,2,:)-XYZ_Atom_Rot(:,1,:);
Vec_CO = squeeze(Vec_CO);
Vec_CO = bsxfun(@rdivide,Vec_CO,sqrt(sum(abs(Vec_CO).^2,2))); % normaliz CO vectors

Vec_CN = XYZ_Atom_Rot(:,3,:)-XYZ_Atom_Rot(:,1,:);
Vec_CN= squeeze(Vec_CN);
Vec_CN = bsxfun(@rdivide,Vec_CN,sqrt(sum(abs(Vec_CN).^2,2))); % normaliz CN vectors

AmideICenter = squeeze(XYZ_Atom_Rot(:,1,:)) + Vec_CO.*0.665 + Vec_CN.*0.256; % center of amide I mode, ref from Jenny's mathematica code

% Define Lab frame coordinate of each mode
Z_Sim = Vec_CO;
X_Sim = cross(Vec_CN,Z_Sim,Dimenesion);
X_Sim = bsxfun(@rdivide,X_Sim,sqrt(sum(abs(X_Sim).^2,Dimenesion))); % normalize
Y_Sim = cross(Z_Sim,X_Sim,Dimenesion);
Y_Sim = bsxfun(@rdivide,Y_Sim,sqrt(sum(abs(Y_Sim).^2,2))); % normalize

XYZ_Sim = [X_Sim(:); Y_Sim(:); Z_Sim(:)]';
XYZ_Sim = reshape(XYZ_Sim,Num_Modes,3,[]); 

%% Calculate the transition dipoles (mu) and Raman tensors (alpha) in the Lab frame

% Transition dipole Angle
mu_angle = -27.5*pi/180; % cunter clockwise roation about x axis from Jenny's paper (10.1021/jp408064b)
% mu_angle = -20*pi/180; % from D. Strasfeld's paper (10.1021/jp9072203)
% mu_angle = -25*pi/180; % from Cho's paper
% mu_angle = -10*pi/180; % follow Skinner's 2014 paper, doi:10.1063/1.4882059

mu_vec0 = 16.1016*[0,0,-1]';       % this scale factor was calculated from DFT simulation, check the 'NMA' note 
mu_Mol  = (Rx(mu_angle)*mu_vec0)'; % 3 by 1

% Raman tensor
alpha_angle = -34*pi/180;    % counter clockwise roation about x axis from Jenny's paper (10.1021/jp408064b)
alpha0 = 5*[ 0.05, 0  , 0;   
             0   , 0.2, 0; 
             0   , 0  , 1 ];

alpha_Mol = Rx(alpha_angle) * alpha0 * Rx(alpha_angle)';

mu_Sim    = zeros(Num_Modes,3);
alpha_Sim = zeros(Num_Modes,3,3);

for ii=1:Num_Modes
    mu_Sim(ii,:,:)    = squeeze(XYZ_Sim(ii,:,:))*mu_Mol';
    alpha_Sim(ii,:,:) = squeeze(XYZ_Sim(ii,:,:))*alpha_Mol'*squeeze(XYZ_Sim(ii,:,:))';
end    

% take vectorize Alpha
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

%% Define Mode frequency and anharmonicity

AmideIFreq = ones(Num_Modes,1)* NLFreq;

if ~ischar(L_Index)
    AmideIFreq(L_Index) = LFreq.*ones(size(L_Index));
end

% anharmonicity = 12
AmideIAnharm = ones(Num_Modes,1).*Anharm;

%% Output Structure
Output.center         = AmideICenter;
Output.freq           = AmideIFreq;
Output.anharm         = AmideIAnharm;
Output.mu             = mu_Sim;
Output.alpha          = alpha; % raman tensor vector form [N x 9]
Output.alpha_matrix   = alpha_Sim;
Output.AtomSerNo      = AmideIAtomSerNo;
Output.Num_Modes      = Num_Modes;
Output.XYZ            = XYZ_Rot;
Output.FilesName      = FilesName;
Output.mu_angle       = mu_angle;
Output.alpha_angle    = alpha_angle;
Output.Mol_Frame_Rot  = Mol_Frame_Rot;
Output.Mol_Frame_Orig = Mol_Frame_Orig;


