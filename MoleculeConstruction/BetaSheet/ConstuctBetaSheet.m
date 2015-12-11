function XYZ = ConstuctBetaSheet(Num_Residue,Num_Strand,TransV,RotV)
%% BetaSheet(sheetSize,E,theta,anhar,du,graph,hpar)
% 
% sheetsize = [# of strands, # of residue per strand]
%  
%   Given sheet size, this script will generate an ideal beta sheet 

% ------- Version log -----------------------------------------------------
% 
% Ver. 1.0  150916  Modified from Betasheet in 2DIR_Betasheet project
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2015

% ---[Debug]--------------------------------------------------------------
% clear all
% Num_Residue= 10;
% Num_Strand = 10;
% ------------------------------------------------------------------------

% initail setup from Lauren's code
% rotate = [0,0,0];
% hpar = [0 0 4.75 rotate(1)*(pi/180) rotate(2)*(pi/180) rotate(3)*(pi/180)];
hpar = [TransV,RotV.*pi./180];
% hpar = [0 0 4.00 rotate(1)*(pi/180) rotate(2)*(pi/180) rotate(3)*(pi/180)];
% E = ones(3*4);
% theta = 20/180*pi;
% du = 3.144;
% anhar = 14;
% graph = 'n';


%%Pull coordinates from short ideal parallel beta-sheet
fid = fopen('BETA-3.ENT');
Template = textscan(fid,'%s%d%s%s%d%f%f%f');
fclose(fid);

AtomName = Template{3};
XYZ_0 = cell2mat(Template(:,6:8));
NumAtoms=length(AtomName);

%% Pull coordinates for each atom of the ideal beta-sheet

N_0  = [];
C_0  = [];
O_0  = [];
CA_0 = [];

for i = 1:NumAtoms
    switch AtomName{i}
        case 'C'
            C_0  = [C_0 ; XYZ_0(i,:)];
        case 'O'
            O_0  = [O_0 ; XYZ_0(i,:)];
        case 'N'
            N_0  = [N_0 ; XYZ_0(i,:)];
        case 'CA'
            CA_0 = [CA_0; XYZ_0(i,:)];
    end
end

%% Calculate average distances and vectors

%Between residues i and i+2 along the peptide backbone
Chain1 = N_0(3,:)-N_0(1,:);
Chain2 = C_0(3,:)-C_0(1,:);
Chain = 0.5*(Chain1+Chain2);
% ChainMag = norm(Chain); % chainMag = sqrt(diag(chain*chain'));
% ChainUnit = Chain/ChainMag;

%Between C and O along carbonyl bond
Bond1 = O_0(1,:)-C_0(1,:);
Bond2 = C_0(2,:)-O_0(2,:);
Bond = 0.5*(Bond1+Bond2);
% BondMag = norm(Bond); % BondMag = sqrt(diag(Bond*Bond'));
% BondUnit = Bond/BondMag;

%% Determine positions of atoms for first strand

C_XYZ=zeros(Num_Residue,3); 
O_XYZ=zeros(Num_Residue,3); 
N_XYZ=zeros(Num_Residue,3); 
CA_XYZ=zeros(Num_Residue,3);

% chain=chainUnit*6.5;  %why round???
for i=1:Num_Residue
    if mod(i,2)==1          %Odd residues
        C_XYZ (i,:) = C_0 (1,:)+floor(i/2)*Chain;
        O_XYZ (i,:) = O_0 (1,:)+floor(i/2)*Chain;
        N_XYZ (i,:) = N_0 (1,:)+floor(i/2)*Chain;
        CA_XYZ(i,:) = CA_0(1,:)+floor(i/2)*Chain;
    else                    %Even residues
        C_XYZ (i,:) = C_0 (2,:)+(i/2-1)*Chain;
        O_XYZ (i,:) = O_0 (2,:)+(i/2-1)*Chain;
        N_XYZ (i,:) = N_0 (2,:)+(i/2-1)*Chain;
        CA_XYZ(i,:) = CA_0(2,:)+(i/2-1)*Chain;
    end
end

%% Find center of mass (CoMass) of strand

% AW => AtomicWeight
AW_C = 12.011;
AW_N = 14.007;
AW_O = 15.999;

% MW => MassWeighted Coordinate
MW_C  =  C_XYZ*AW_C;
MW_N  =  N_XYZ*AW_N;
MW_O  =  O_XYZ*AW_O;
MW_CA = CA_XYZ*AW_C;

Strand_Weight = Num_Residue*(2*AW_C + AW_O + AW_N);

CoMass = sum((MW_C + MW_N + MW_O + MW_CA),1)/Strand_Weight;


%% Add in proper rotation and angles of residues along strand

gpar=[1.2165 3.8888 -2.4406];   % What is this?
nma_md=[0 0 0; 0 0 4.75];       % Z distance? 

nma_es(1,:)=CoMass;             % NMA estimation?
nma_es(2,:)=CoMass+gpar;

for nn=1:10
    if nn==1
        fpar(1:3)=[14.5004 -9.4375 -6.3790];    % Translation parameter
        fpar(4:6)=[-4.0559 -0.5700 0.6973];     % Rotation parameter in radius
    else
        fpar=out;
    end
    options=optimset('MaxFunEvals',10000);
    [xx1,fval,exitflag]=fminsearch('rotate_nma',fpar,options,nma_es,nma_md);
    out=xx1;
end
fpar=out;

Rx = [1 0 0; 0 cos(fpar(4)) -sin(fpar(4)); 0 sin(fpar(4)) cos(fpar(4))];
Ry = [cos(fpar(5)) 0 sin(fpar(5)); 0 1 0; -sin(fpar(5)) 0 cos(fpar(5))];
Rz = [cos(fpar(6)) -sin(fpar(6)) 0; sin(fpar(6)) cos(fpar(6)) 0; 0 0 1];

% Create dummy variables and apply rotation and translation to first strand
ca2=zeros(size(CA_XYZ)); n2=zeros(size(N_XYZ)); c2=zeros(size(C_XYZ)); o2=zeros(size(O_XYZ));
for i=1:length(CA_XYZ)
    ca2(i,:) = (Rx*Ry*Rz*(CA_XYZ(i,:)' + fpar(1:3)'))';
    n2(i,:) = (Rx*Ry*Rz*(N_XYZ(i,:)' + fpar(1:3)'))';
    c2(i,:) = (Rx*Ry*Rz*(C_XYZ(i,:)' + fpar(1:3)'))';
    o2(i,:) = (Rx*Ry*Rz*(O_XYZ(i,:)' + fpar(1:3)'))';
end

clear coord carbonA carbon nitrogen oxygen
CA_XYZ=ca2; N_XYZ=n2; C_XYZ=c2; O_XYZ=o2;
clear ca2 n2 c2 o2

%% Group all atoms in first strand together and find COM of each residue

for i=1:Num_Residue
    coord1((i-1)*4+1,:)=CA_XYZ(i,:);
    coord1((i-1)*4+2,:)=N_XYZ(i,:);
    coord1((i-1)*4+3,:)=C_XYZ(i,:);    
    coord1((i-1)*4+4,:)=O_XYZ(i,:);
    temp1=CA_XYZ(i,:)*12.011;
    temp2=N_XYZ(i,:)*14.007;
    temp3=C_XYZ(i,:)*12.011;
    temp4=O_XYZ(i,:)*15.999;
    COM(i,:)=(temp1+temp2+temp3+temp4)/(12.011+14.007+12.011+15.999);
end

clear carbonA nitrogen carbon oxygen temp1 temp2 temp3 temp4

%% Translate first strand to create the full beta-sheet

T=[hpar(1) hpar(2) hpar(3)];
Rx = [1 0 0; 0 cos(hpar(4)) -sin(hpar(4)); 0 sin(hpar(4)) cos(hpar(4))];
Ry = [cos(hpar(5)) 0 sin(hpar(5)); 0 1 0; -sin(hpar(5)) 0 cos(hpar(5))];
Rz = [cos(hpar(6)) -sin(hpar(6)) 0; sin(hpar(6)) cos(hpar(6)) 0; 0 0 1];

for i=1:Num_Strand
    for nn=1:length(coord1)
        if i==1
            coord2((i-1)*length(coord1)+nn,:)=coord1(nn,:);
        else
            coord2((i-1)*length(coord1)+nn,:)=(coord2((i-2)*length(coord1)+nn,:)+T)*Rx*Ry*Rz;
        end
    end
end

CA_XYZ=coord2(1:4:length(coord2),:);
N_XYZ=coord2(2:4:length(coord2),:);
C_XYZ=coord2(3:4:length(coord2),:);
O_XYZ=coord2(4:4:length(coord2),:);

clear coord1 coord2

%% Rotate the whole molecule so strand align with X axis 
Strand_Axis = C_XYZ(end,:)-C_XYZ(1,:);
Strand_Axis = Strand_Axis./norm(Strand_Axis);
R_Angle = acos(dot(Strand_Axis,[1,0,0]));
RM = R1_ZYZ_0(R_Angle,0,0);

C_XYZ = (RM*C_XYZ')';
O_XYZ = (RM*O_XYZ')';
N_XYZ = (RM*N_XYZ')';

%% Formating output coordinate
% XYZ = reshape([C_XYZ;O_XYZ;N_XYZ],[],3,3);
XYZ(:,:,1) = C_XYZ;
XYZ(:,:,2) = O_XYZ;
XYZ(:,:,3) = N_XYZ;


% StrucInfo = GetAmideI_CON_XYZ_version(XYZ);
% 
% 
% % Plot FTIR
% 
% PlotFTIR(LocModeNum,mu01Ex,Freq01Ex)
