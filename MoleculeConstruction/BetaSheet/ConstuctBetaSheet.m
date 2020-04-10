function S_BSheet = ConstuctBetaSheet(GUI_Inputs)
%% ConstuctBetaSheet(SheetType,N_Residue,N_Strand,TransV,TwistV)
% Given inputs, ConstuctBetaSheet build idea betasheet so that the strand
% aligned with X axis and the C=O axis alignd with Z axis while its center
% of mass localte at coordinate origin. The betasheet is place in this way
% so the following mode determining code (i.e. GetAmideI.m) could properly
% tranlate/rotate the molecule and its corresponding transition dipoles/
% Raman tensors in lab frame.
%   
% Inputs:
% -------------------------------------------------------------------------
%   sheet type: Anti-parallel(Anit) or Parallele batesheet(Para)
%   # of residue per strand
%   # of strands, , 
%   Strand translational vector (Move along X,Y,Z axis), 
%   Strans twist angles (Twist around X,Y,Z axis)

% Outputs:
%   Num_Atoms: number of atoms that determine how manay loops the mode
%              determine code (i.e. GetAmideI.m)is going to run.
%   XYZ:       XYZ position of each atom.
%   AtomName:  pdb style Atom name, so the mode determining code can tell
%              whay kind of atom is it. (e.g. C, CA, etc...)
% -------------------------------------------------------------------------

% ------- Version log -----------------------------------------------------
% 
% Ver. 1.0  160130  Modified from ConstuctBetaSheet_old
% 
% -------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2016

% ---[Debug]--------------------------------------------------------------
% GUI_Inputs.Debug = 'Debug';
% ------------------------------------------------------------------------

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultSheetType = 'Anti-Parallel';
defaultN_Residue  = 3;
defaultN_Strand   = 2;
defaultTrans_X    = 0;
defaultTrans_Y    = 0;
defaultTrans_Z    = 0;
defaultTwist_X    = 0;
defaultTwist_Y    = 0;
defaultTwist_Z    = 0;
defaultPhi_D      = 0;
defaultPsi_D      = 0;
defaultTheta_D    = 0;
defaultNLFreq     = 1650;
defaultAnharm     = 10;
defaultLFreq      = 1630;
defaultL_Index    = 'None';

% add Optional inputs / Parameters
addOptional(INPUT,'SheetType', defaultSheetType );
addOptional(INPUT,'N_Residue', defaultN_Residue );
addOptional(INPUT,'N_Strand' , defaultN_Strand  );
addOptional(INPUT,'Trans_X'  , defaultTrans_X   );
addOptional(INPUT,'Trans_Y'  , defaultTrans_Y   );
addOptional(INPUT,'Trans_Z'  , defaultTrans_Z   );
addOptional(INPUT,'Twist_X'  , defaultTwist_X   );
addOptional(INPUT,'Twist_Y'  , defaultTwist_Y   );
addOptional(INPUT,'Twist_Z'  , defaultTwist_Z   );
addOptional(INPUT,'Phi_D'    , defaultPhi_D     );
addOptional(INPUT,'Psi_D'    , defaultPsi_D     );
addOptional(INPUT,'Theta_D'  , defaultTheta_D   );
addOptional(INPUT,'NLFreq'   , defaultNLFreq    );
addOptional(INPUT,'Anharm'   , defaultAnharm    );
addOptional(INPUT,'LFreq'    , defaultLFreq     );
addOptional(INPUT,'L_Index'  , defaultL_Index   );

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
SheetType = INPUT.Results.SheetType;
N_Residue = INPUT.Results.N_Residue;
N_Strand  = INPUT.Results.N_Strand;
Trans_X   = INPUT.Results.Trans_X;
Trans_Y   = INPUT.Results.Trans_Y;
Trans_Z   = INPUT.Results.Trans_Z;
Twist_X   = INPUT.Results.Twist_X;
Twist_Y   = INPUT.Results.Twist_Y;
Twist_Z   = INPUT.Results.Twist_Z;
Phi_D     = INPUT.Results.Phi_D;
Psi_D     = INPUT.Results.Psi_D;
Theta_D   = INPUT.Results.Theta_D;
NLFreq    = INPUT.Results.NLFreq;
Anharm    = INPUT.Results.Anharm;
LFreq     = INPUT.Results.LFreq;
L_Index   = INPUT.Results.L_Index;

%% Post process of inputs
TransV = [Trans_X,Trans_Y,Trans_Z];
TwistV = [Twist_X,Twist_Y,Twist_Z];

%% Pull coordinates from short ideal parallel beta-sheet
fid = fopen('BETA-3.ENT');
Template = textscan(fid,'%s%d%s%s%d%f%f%f');
fclose(fid);

AtomName_0  = Template{3};
XYZ_0       = cell2mat(Template(:,6:8));
% Num_Atoms_0 = length(AtomName_0);

%% Rotate the template so that the Chain aixs aligned with X and Bond axis aligned with Z
% Define Chain Vector (X); and Bond Vector (X)
%                 O(6)
%                 ||
% C(1)-N(3)-CA(4)-C(5)-N(7)-CA(8)-C(9)-N(11)-CA(12)
% ||                              ||
% O(2)                            O(10)

%Between residues 1 and 3 along the peptide backbone
Chain1 = XYZ_0(11,:) - XYZ_0(3,:);
Chain2 = XYZ_0( 9,:) - XYZ_0(1,:);
ChainV = 0.5*(Chain1+Chain2);
ChainL = norm(ChainV);

%Between C and O along carbonyl bond
Bond1 = XYZ_0(2,:) - XYZ_0(1,:);
Bond2 = XYZ_0(5,:) - XYZ_0(6,:);
BondV = 0.5*(Bond1+Bond2);

X_Axis = ChainV./norm(ChainV);
Z_Axis = BondV./norm(BondV);
Y_Axis = cross(Z_Axis,X_Axis);
Y_Axis = Y_Axis./norm(Y_Axis);

Target_frame   = [ 1,0,0; 0,1,0; 0,0,1]';
Original_frame = [X_Axis;Y_Axis;Z_Axis]';
RM = Euler_Rot(Target_frame,Original_frame);

XYZ_0_LabFrame = (RM * XYZ_0')'; 
% Reshape so atom split into layers [4X3X3]
% dim1: Residue1, Residue2,Residue3
% dim2: X,Y,Z
% dim3: C,O,N,CA
XYZ_0_LabFrame = permute(reshape(XYZ_0_LabFrame,4,[],3),[2,3,1]);

%% Extend residues along the 1st strand
XYZ_1strand = zeros(N_Residue,3,4);
AtomName_1strand = cell(N_Residue,4);
ChainV_LabFrame = ChainL.*[1,0,0]';

for i = 1:N_Residue
    if mod(i,2)==1
        XYZ_1strand(i,:,:) = bsxfun(@plus,squeeze(XYZ_0_LabFrame(1,:,:)),  floor(i/2).*ChainV_LabFrame);
    else 
        XYZ_1strand(i,:,:) = bsxfun(@plus,squeeze(XYZ_0_LabFrame(2,:,:)),floor(i/2-1).*ChainV_LabFrame);
    end
    
    AtomName_1strand(i,:) = AtomName_0(1:4);
end

AtomName_1strand = reshape(AtomName_1strand',[],1);

%% substitute the last N and CA to O and H on C terminus for the first strand
% add x at the end of the atom name so it does not confuse the GetAmideI.m
Ind_H_1strand = length(AtomName_1strand);
AtomName_1strand{Ind_H_1strand-3} = 'Cx';
AtomName_1strand{Ind_H_1strand-2} = 'Ox';
AtomName_1strand{Ind_H_1strand-1} = 'Ox';
AtomName_1strand{Ind_H_1strand}   = 'Hx';

%% Define flipping center for the fliping operation in creating antiparallel betasheet
if mod(N_Residue,2)
    % N_Residue is odd, but number of amide mode is even
    Ind_Center_AmideI = (N_Residue-1)/2;
    
    COF_O1_XYZ = XYZ_1strand(Ind_Center_AmideI  ,:,2);
    COF_O2_XYZ = XYZ_1strand(Ind_Center_AmideI+1,:,2);
    COF_N1_XYZ = XYZ_1strand(Ind_Center_AmideI  ,:,3);
    COF_N2_XYZ = XYZ_1strand(Ind_Center_AmideI+1,:,3);
    COF_V = (COF_O1_XYZ + COF_O2_XYZ + COF_N1_XYZ + COF_N2_XYZ)./4;
    COF   = reshape(COF_V,1,3,1);

else
    % N_Residue is even, but number of amide mode is odd
    Ind_Center_AmideI = (N_Residue)/2;
    
    COF_C_XYZ = XYZ_1strand(Ind_Center_AmideI,:,1);
    COF_N_XYZ = XYZ_1strand(Ind_Center_AmideI,:,3);
    COF_V = (COF_C_XYZ + COF_N_XYZ)./2;
    COF   = reshape(COF_V,1,3,1);
end

% move the 1st strand's COF to origin of lab fram
XYZ_1strand_COF = bsxfun(@minus,XYZ_1strand,COF);

%% reshape to C,O,N,Ca going along 1st dimension and XYZ along 2nd
% [Cx ,Cy ,Cz ]
% [Ox ,Oy ,Oz ]
% [Nx ,Ny ,Nz ]
% [CAx,CAy,CAz]
% ...
% [Cx ,Cy ,Cz ]
% [Ox ,Oy ,Oz ]
% [Ox ,Oy ,Oz ]
% [Hx ,Hy ,Hz ]

XYZ_1strand_COF = reshape(permute(XYZ_1strand_COF,[3,1,2]),[],3);
N_Atom_1strand = size(XYZ_1strand_COF,1);

%% Generating the rest of strands
XYZ = zeros(N_Atom_1strand,N_Strand,3);
XYZ(:,1,:) = XYZ_1strand_COF;
AtomName = cell(N_Atom_1strand,N_Strand);
AtomName(:,1)= AtomName_1strand;

if strcmp(SheetType,'Anti-Parallel')
    APB_Flag = 1;
    SheetTypeString = 'APB';
else
    APB_Flag = 0;
    SheetTypeString = 'PB';
end

TwistV = TwistV./180.*pi; % unit raidus
for j = 2:N_Strand
    
    if APB_Flag
        if mod(N_Residue,2)
        % Rotate strand around Y axis for APB with odd number of strand
            XYZ_j_tmp1 = (Ry((j-1)*pi) * XYZ_1strand_COF')';
        else
        % flip X to -X of strand for APB with even number of strand
            Flip_V = [(-1)^(j-1),1,1];
            XYZ_j_tmp1 = bsxfun(@times,XYZ_1strand_COF,Flip_V);
        end
        
        XYZ_j_flip = XYZ_j_tmp1;
        AtomName_1strand_j = AtomName_1strand;
    else 
        % for PB doing notheing here
        XYZ_j_flip = XYZ_1strand_COF; 
        AtomName_1strand_j = AtomName_1strand;
    end
    
    % twist strand
    %TwistM = (Rx(TwistV(1))*Ry(TwistV(2))*Rz(TwistV(3)) )^(j-1);
    TwistM = ( R1_ZYZ_0(TwistV(1),TwistV(2),TwistV(3)) )^(j-1);
    XYZ_j_flip_tw = (TwistM * XYZ_j_flip')';
    
    % move strand
    XYZ(:,j,:) = bsxfun(@plus, XYZ_j_flip_tw, (j-1)*TransV);  
    
    % Save Atome name for each strand
    AtomName(:,j)= AtomName_1strand_j;
end

% reshape
AtomName = reshape(AtomName,[],1);
XYZ      = reshape(XYZ,[],3);

%% Retrieve O and H index of C terminus for molecule plotting
Num_Atoms = length(AtomName);
Atoms_Array = 1:Num_Atoms;
Ind_H = Atoms_Array(strcmp(AtomName,'H'));
Ind_O = Ind_H -1;

%% Initiate a StrutureData class container to fill in information
S_BSheet = StructureData;
S_BSheet.XYZ        = XYZ;
S_BSheet.AtomName   = AtomName;
S_BSheet.Extra.ChainID = {'     '}; % output dummy space array to skip Chain recognition 

% Find all the amide modes and their default values
S_BSheet = SD_GetAmideI(S_BSheet);

%% Rotate the molecule
RotV_D   = [Phi_D,Psi_D,Theta_D];
RotV     = RotV_D./180.*pi;
R        = R1_ZYZ_0(RotV(1),RotV(2),RotV(3));
S_BSheet = SD_Rot(S_BSheet,R);

%% Collecting Extra info
Extra = S_BSheet.Extra;
Extra.N_Residue         = N_Residue;
Extra.N_Strand          = N_Strand;
Extra.N_Mode_per_Starnd = N_Residue-1;
Extra.Ind_H             = Ind_H;
Extra.Ind_O             = Ind_O;
Extra.TransV            = TransV;
Extra.TwistV            = TwistV;
Extra.RotV_D            = RotV_D;
S_BSheet.Extra = Extra;

%% Modify Mode frequency, anharmonicity, and Labling according to the GUI inputs
LocFreq   = ones(S_BSheet.Nmodes,1).*NLFreq;
LocAnharm = ones(S_BSheet.Nmodes,1).*Anharm;

if ~ischar(L_Index)
    LocFreq(L_Index) = LFreq;
end

S_BSheet.LocAnharm  = LocAnharm;
S_BSheet.LocFreq    = LocFreq;

%% Formating other outputs 
S_BSheet.FilesName  = [SheetTypeString,'-R',num2str(N_Residue),'S',num2str(N_Strand)];
S_BSheet.GUI_Inputs = GUI_Inputs;
S_BSheet.hPlotFunc  = @Plot_Betasheet_AmideI; % Add figure drawing function and GUI inputs
S_BSheet            = SD_1ExH(S_BSheet); % Calculate One Exciton Hamiltoian
