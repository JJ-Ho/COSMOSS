function Output = ConstuctBetaSheet(SheetType,N_Residue,N_Strand,TransV,TwistV)
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
% clear all
% SheetType = 'Anti';
% N_Residue= 7;
% N_Strand = 1;
% TransV = [0,0,4.75];
% TwistV = [0,0,4];
% ------------------------------------------------------------------------

%% Pull coordinates from short ideal parallel beta-sheet

% hpar = [TransV,TwistV.*pi./180];

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
N_Atom_per_strand = length(AtomName_1strand);
%% substitute the last N and CA to O and H on C terminus
Ind_H_1strand = length(AtomName_1strand);
Ind_O_1strand = length(AtomName_1strand)-1;
AtomName_1strand{Ind_O_1strand} = 'O';
AtomName_1strand{Ind_H_1strand} = 'H';

Ind_H = Ind_H_1strand:N_Atom_per_strand:N_Strand*N_Atom_per_strand;
Ind_O = Ind_O_1strand:N_Atom_per_strand:N_Strand*N_Atom_per_strand-1;
%% Define flipping center for the fliping operation in creating antiparallel betasheet
if mod(N_Residue,2)
    % N_Residue is odd, but number of amide mode is even
    Ind_Center_AmideI = (N_Residue-1)/2+1;
else
    % N_Residue is even, but number of amide mode is odd
    Ind_Center_AmideI = (N_Residue)/2;
end

COF_C_XYZ = XYZ_1strand(Ind_Center_AmideI,:,1);
COF_N_XYZ = XYZ_1strand(Ind_Center_AmideI,:,3);
COF_V = (COF_C_XYZ + COF_N_XYZ)./2;
COF   = reshape(COF_V,1,3,1);
    
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

if strcmp(SheetType,'Anti')
    Flip_Sign = -1;
    SheetTypeString = 'APB';
else
    Flip_Sign = 1;
    SheetTypeString = 'PB';
end

TwistV = TwistV/180*pi; % unit raidus
for j = 2:N_Strand
    
    % flip strand if Anti-parallel
    Flip_V = [Flip_Sign^(j-1),1,1];
    XYZ_j_flip = bsxfun(@times,XYZ_1strand_COF,Flip_V);
    
    % twist strand
    TwistM = (Rx(TwistV(1))*Ry(TwistV(2))*Rz(TwistV(3)) )^(j-1);
    XYZ_j_flip_tw = (TwistM * XYZ_j_flip')';
    
    % move strand
    XYZ(:,j,:) = bsxfun(@plus, XYZ_j_flip_tw, (j-1)*TransV);  
    
    % Save Atome name for each strand
    AtomName(:,j)= AtomName_1strand;
end

% reshape
AtomName = reshape(AtomName,[],1);
XYZ      = reshape(XYZ,[],3);

%% Formating output coordinate
Output.Num_Atoms = length(AtomName);
Output.XYZ       = XYZ;
Output.AtomName  = AtomName;
Output.FilesName = [SheetTypeString,'-R',num2str(N_Residue),'S',num2str(N_Strand)];
Output.Ind_H     = Ind_H;
Output.Ind_O     = Ind_O;
