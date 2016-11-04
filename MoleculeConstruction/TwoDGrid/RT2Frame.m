function Output = RT2Frame(S_Info)

%% Process Inputs
XYZ      = S_Info.XYZ;
TDV      = S_Info.TDV;
Raman    = S_Info.Raman;
Mode_Num = S_Info.Mode_Num;


Mol_Frame  = S_Info.Mol_Frame;
Center_Ind = Mol_Frame.Center_Ind;
Z_i_Ind    = Mol_Frame.Z_i_Ind;
Z_f_Ind    = Mol_Frame.Z_f_Ind;
XY_i_Ind   = Mol_Frame.XY_i_Ind;
XY_f_Ind   = Mol_Frame.XY_f_Ind;
Frame_Type = Mol_Frame.Frame_Type;

%% Get Center and Raotation matrix
% Get Center XYZ
C = sum(XYZ(Center_Ind,:),1)/length(Center_Ind);
XYZ_T = bsxfun(@minus,XYZ,C);

% Get Ratational Matrix to rotate the molecule to defined frame 
Vec_Z  = XYZ(Z_f_Ind,:) - XYZ(Z_i_Ind,:);
Z      = Vec_Z/norm(Vec_Z);
Vec_XZ = XYZ(XY_f_Ind,:) - XYZ(XY_i_Ind,:);
Vec_XZ = Vec_XZ/norm(Vec_XZ);
Y      = cross(Z,Vec_XZ);
Y      = Y/norm(Y);
X      = cross(Y,Z);
X      = X/norm(X);

Mol_Frame      =zeros(3,3);
Mol_Frame(:,1) = X;
Mol_Frame(:,2) = Y;
Mol_Frame(:,3) = Z;
New_Frame      = eye(3);
Rot = Euler_Rot(New_Frame,Mol_Frame);

% Rotate molecule definition from X-Z plane to Y-Z plane
if strcmp(Frame_Type,'YZ')
    R_XZ2YZ = R1_ZYZ_0(-pi/2,0,0);
    Rot = R_XZ2YZ*Rot; 
end

XYZ_T_R = (Rot*XYZ_T')';

%% Rotate the TDV and Raman tensor accordingly
% Rotated Transition Dipole Vector
TDV_Rot = (Rot*TDV')';

% Rotated Raman Tensor
Raman_Rot = zeros(size(Raman));
for i = 1:Mode_Num
    Raman_Rot(:,:,i) = Rot'*squeeze(Raman(:,:,i))*Rot;
end

%% Output
Output.C_MF         = C;
Output.Rot2MF       = Rot;

Output.MF       = S_Info;
Output.MF.XYZ   = XYZ_T_R;
Output.MF.TDV   = TDV_Rot;
Output.MF.Raman = Raman_Rot;

Output.Orig = S_Info;