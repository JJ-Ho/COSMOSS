function P = ReadG09Input(FName,RR,varargin)
% 
% Given a input file name "FName", this script can parse XYZ, mu and alpha
% variables and generate parsed input structure "P"
% 

% ------- Version log ----------------------------------------------------
% 
% V1.1  150426  Copy from 2DSFG_Ester project
%               Update debug part
% 
% V1.0  140903  Modified from ReadDNAInput v1.6
%               Add Conn(3,13) = 1 for connection S-C bond
% 
% ------------------------------------------------------------------------ 
% Copyright Jia-Jung Ho, 2013

% ----[Debug] ------------------------------
% clear all
% close all
% clc
% FName = '140903_MMB.txt';
% RR = [0,0,0];
% varargin = {'MolFrame','XZ',...
%             };
% ------------------------------------------

%% Input part

Inputs = inputParser;

% Default values
defaultMolFrame   = 'XZ'; 
expectedMolFrame   = {'XZ','YZ'};

% Add optional inputs to inputparser object
addParamValue(Inputs,'MolFrame',defaultMolFrame,...
                 @(x) any(validatestring(x,expectedMolFrame)));

parse(Inputs,varargin{:});

% Reassign Variable names
MolFrame   = Inputs.Results.MolFrame;

%% Start to read inputs
fid = fopen(FName);

%% Predefined Variables
Atom_Num = cell2mat(textscan(fid,'[Atom_Num] %f',1,'commentStyle','%'));
Mode_Num = cell2mat(textscan(fid,'[Mode_Num] %d',1,'commentStyle','%'));

%% XYZ part
Atom = textscan(fid,'[Atom] %s %f %f %f',Atom_Num,'CollectOutput',1,'commentStyle','%');

Atom_Name = Atom{1}; % Atom symbel
XYZ_Orig = Atom{2};  % Unrotated molecule xyz. Already a array since 'CollectOutput'

% Connectivity
[n,m] = ndgrid(1:Atom_Num);
T1=XYZ_Orig(m(:),:)-XYZ_Orig(n(:),:);
T2=sqrt(sum(T1.^2,2));
Distance_matrix=reshape(T2,Atom_Num,Atom_Num);
lower=tril(Distance_matrix,-1);
[a,b]=find(lower<1.6 & lower>0);
C_index=[a,b];
Conn=false(Atom_Num);
Conn(C_index(:,1)+(C_index(:,2)-1)*Atom_Num)=true;

Conn(13,3) = 1; % Add S-C connection for Ester

Conn=Conn|Conn';

%% Molecular orientation part + Mu part
% [Orientation]  Center Z_I Z_F XZ_I XZ_F
Orientation = textscan(fid,'[Orientation] %f %f %f %f %f',1,'commentStyle','%','CollectOutput',1);
Orientation = cell2mat(Orientation);

% Considering the case that we want to rotate alone mu instead of molecule 
% in this case the [orientation] will be 0,0,0,0,0
if any(Orientation)    

    Center = XYZ_Orig(Orientation(1),:);
    Vec_Z  = XYZ_Orig(Orientation(3),:) - XYZ_Orig(Orientation(2),:);
    Z = Vec_Z/norm(Vec_Z);
    Vec_XZ = XYZ_Orig(Orientation(5),:) - XYZ_Orig(Orientation(4),:);
    Vec_XZ = Vec_XZ/norm(Vec_XZ);

    XYZ_T = bsxfun(@minus,XYZ_Orig,Center);

    Y = cross(Z,Vec_XZ);
    Y = Y/norm(Y);
    X = cross(Y,Z);
    X = X/norm(X);

    Mol_Frame=zeros(3,3);
    Mol_Frame(:,1) = X;
    Mol_Frame(:,2) = Y;
    Mol_Frame(:,3) = Z;
    New_Frame = eye(3);
    ERot = Euler_Rot(New_Frame,Mol_Frame);

    % Rotate molecule definition from X-Z plane to Y-Z plane
    if strcmp(MolFrame,'YZ')
        R_XZ2YZ = R1_ZYZ_0(-pi/2,0,0);
        ERot = R_XZ2YZ*ERot; 
    end
    
    R_Asigned = R1_ZYZ_0(RR(1),RR(2),RR(3));
    R_Total   = R_Asigned*ERot;
    
    XYZ_T_R = (R_Total*XYZ_T')';
    
    % Mu Vector Part
    % Unrotated Transition Dipole Vector (mu)
    TDV_Orig = textscan(fid,'[TDV] %s %f %f %f',Mode_Num,'CollectOutput',1,'commentStyle','%');
    TDV_Orig = TDV_Orig{2};
    % Rotated Transition Dipole Vector
    TDV_Rot = (R_Total*TDV_Orig')';
    
else
    % if all element of orientation are zeros => rotate transition dipole 
    ring_ind = 1:6;
    Center = sum(XYZ_Orig(ring_ind,:),1)./length(ring_ind);
    XYZ_T = bsxfun(@minus,XYZ_Orig,Center);
    
    % Read-in Vector sequence: [Vec_Seq] VecZ VecXZ
    Vec_Seq = textscan(fid,'[Vec_Seq] %d %d',1,'CollectOutput',1,'commentStyle','%');
    Vec_Seq = cell2mat(Vec_Seq);
    
    % Unrotated Transition Dipole Vector (mu)
    TDV_Orig = textscan(fid,'[TDV] %s %f %f %f',Mode_Num,'CollectOutput',1,'commentStyle','%');
    TDV_Orig = TDV_Orig{2};

    Vec_Z  = TDV_Orig(Vec_Seq(1),:);
    Z = Vec_Z/norm(Vec_Z);
    Vec_XZ = TDV_Orig(Vec_Seq(2),:);
    Vec_XZ = Vec_XZ/norm(Vec_XZ);
    
    Y = cross(Z,Vec_XZ);
    Y = Y/norm(Y);
    X = cross(Y,Z);
    X = X/norm(X);
    
    Mol_Frame=zeros(3,3);
    Mol_Frame(:,1) = X;
    Mol_Frame(:,2) = Y;
    Mol_Frame(:,3) = Z;
    New_Frame = eye(3);
    ERot = Euler_Rot(New_Frame,Mol_Frame);
    
    % Rotate molecule definition from X-Z plane to Y-Z plane
    if strcmp(MolFrame,'YZ')
        R_XZ2YZ = R1_ZYZ_0(-pi/2,0,0);
        ERot = R_XZ2YZ*ERot; 
    end
    
    R_Asigned = R1_ZYZ_0(Orientation(1),Orientation(2),Orientation(3));
    R_Total   = R_Asigned*ERot;
    
    XYZ_T_R = (R_Total*XYZ_T')';  
    TDV_Rot = (R_Total*TDV_Orig')';
end

%% Alpha Vector Part
% Unrotated Raman Tensor (alpha)
Raman_Orig = textscan(fid,'[Raman] %s %f %f %f %f %f %f %f %f %f',Mode_Num,'CollectOutput',1,'commentStyle','%');
Raman_Matrix = reshape((Raman_Orig{2})',[3,3,Mode_Num]);

% Rotated Raman Tensor
Raman_Rot = zeros(size(Raman_Matrix));
for i = 1:Mode_Num
    Raman_Rot(:,:,i) = R_Total*squeeze(Raman_Matrix(:,:,i))*R_Total';
end

% Vectorize Raman tensor with following convention
% [ Aixx Aixy Aixz ] 
% [ Aiyx Aiyy Aiyz ]
% [ Aizx Aizy Aizz ]
% 
% [Aixx Aiyx Aizx Aixy Aiyy Aizy Aixz Aiyz Aizz]
% 

Raman_Rot_Vec = reshape(Raman_Rot,[9,Mode_Num])';


%% Intensity scaling of IR and Raman

Int_Harm.IR      = textscan(fid,'[Int_Harm]   IR    %s %f',Mode_Num,'commentStyle','%');
Int_Harm.Raman   = textscan(fid,'[Int_Harm]   Raman %s %f',Mode_Num,'commentStyle','%');
Int_AnHarm.IR    = textscan(fid,'[Int_AnHarm] IR    %s %f',Mode_Num,'commentStyle','%');
Int_AnHarm.Raman = textscan(fid,'[Int_AnHarm] Raman %s %f',Mode_Num,'commentStyle','%');


%% Frequency Part
% Toggle of Anharmonic correction
AnharmCorrect   = textscan(fid,'[AnharmCorrect] %s %f',1,'commentStyle','%');
FreqScaleFactor = textscan(fid,'[FreqScaleFactor] %f',1,'commentStyle','%'); 
FreqScaleFactor = cell2mat(FreqScaleFactor);

if strcmp(AnharmCorrect{1},'no')
    Anharm = AnharmCorrect{2};
    Freq_tmp = textscan(fid,'[Freq] %s %f',Mode_Num,'commentStyle','%');
    Fundamental = Freq_tmp{2};
    Overtone = 2*Fundamental-Anharm;
    
    % Create index to calculate combination band
    [I2,I1] = ndgrid(1:Mode_Num,1:Mode_Num);
    Ind_1 = I1(logical(tril(ones(3,3),-1)));
    Ind_2 = I2(logical(tril(ones(3,3),-1)));
    Combination = Fundamental(Ind_1)+Fundamental(Ind_2)-Anharm;
    
elseif strcmp(AnharmCorrect{1},'yes')
    Fundamental = textscan(fid,'[Freq] %s %f',Mode_Num,'commentStyle','%');
    Overtone    = textscan(fid,'[Freq] %s %f',Mode_Num,'commentStyle','%');
    Combination = textscan(fid,'[Freq] %s %f',Mode_Num*(Mode_Num-1)/2,'commentStyle','%');
    
    Fundamental = Fundamental{2};
    Overtone    = Overtone{2};
    Combination = Combination{2};
    
else 
    disp('Please provide correct "[AnharmCorrect]" option [yes|no]...')
end

Freq.Fundamental = FreqScaleFactor*Fundamental;
Freq.Overtone    = FreqScaleFactor*Overtone;
Freq.Combination = FreqScaleFactor*Combination;

fclose(fid);
%% Formate Output

P.FName           = FName;
P.Atom_Num        = Atom_Num;
P.Mode_Num        = Mode_Num;
P.Atom_Name       = Atom_Name;
P.Conn            = Conn;
P.XYZ             = XYZ_T_R;
P.Freq            = Freq;
P.TDV             = TDV_Rot;
P.Raman           = Raman_Rot;
P.RamanVec        = Raman_Rot_Vec;
P.Int_Harm.IR     = Int_Harm.IR{2}   ;
P.Int_Harm.Raman  = Int_Harm.Raman{2};
P.Int_AnHarm.IR   = Int_AnHarm.IR{2}   ;
P.Int_AnHarm.Raman= Int_AnHarm.Raman{2};
P.ERot            = ERot;
P.R_Asigned       = R_Asigned;
P.Trans           = Center;

P.Orig.XYZ              = XYZ_Orig;
P.Orig.TDV              = TDV_Orig;
P.Orig.Raman            = Raman_Orig;
P.Orig.Freq.Fundamental = Fundamental;
P.Orig.Freq.Overtone    = Overtone;
P.Orig.Freq.Combination = Combination;


