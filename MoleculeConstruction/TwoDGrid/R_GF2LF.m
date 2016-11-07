function Output = R_GF2LF(S_Info)
% This function transform structural and mode related information from 
% original coordinate frame (G09 simulation frame) to lab frame so that 
% the rotatable bond aligned with z axis. This function take the atomic
% index from formatted G09 output file or GUI inputs to define the new
% origin and axes in molecule frame. By selection, this function will do
% rotational average around z axis in mol. frame to account for free
% rotatable.

%% Process Inputs
XYZ      = S_Info.XYZ;
TDV      = S_Info.TDV;
Raman    = S_Info.Raman;
Mode_Num = S_Info.Mode_Num;


MF_Info = S_Info.MF_Info;
Center_Ind = MF_Info.Center_Ind;
Z_i_Ind    = MF_Info.Z_i_Ind;
Z_f_Ind    = MF_Info.Z_f_Ind;
XY_i_Ind   = MF_Info.XY_i_Ind;
XY_f_Ind   = MF_Info.XY_f_Ind;
Frame_Type = MF_Info.Frame_Type;
BondAvg    = MF_Info.BondAvg;

LF_Info = S_Info.LF_Info;
LF_Phi   = LF_Info.LF_Phi;
LF_Psi   = LF_Info.LF_Psi;
LF_Theta = LF_Info.LF_Theta;

%% Get Center and Rotation matrix
% Get Center XYZ
C_G09_Frame = sum(XYZ(Center_Ind,:),1)/length(Center_Ind);
XYZ_T = bsxfun(@minus,XYZ,C_G09_Frame);

% Get Ratational Matrix to rotate the molecule to defined frame 
Vec_Z  = XYZ(Z_f_Ind,:) - XYZ(Z_i_Ind,:);
Z      = Vec_Z/norm(Vec_Z);
Vec_XZ = XYZ(XY_f_Ind,:) - XYZ(XY_i_Ind,:);
XZ     = Vec_XZ/norm(Vec_XZ);
Y      = cross(Z,XZ);
Y      = Y/norm(Y);
X      = cross(Y,Z);
X      = X/norm(X);

% Generate rotational matrix from G09 simulation frame to defined molecular
% frame
Mol_Frame      = zeros(3,3);
Mol_Frame(:,1) = X;
Mol_Frame(:,2) = Y;
Mol_Frame(:,3) = Z;
New_Frame      = eye(3);
Rot2MF = Euler_Rot(New_Frame,Mol_Frame);

% Rotate molecule definition from X-Z plane to Y-Z plane
switch Frame_Type
    case 1 %Frame_Type = 'XZ';
    case 2 %Frame_Type = 'YZ';
        R_XZ2YZ = R1_ZYZ_0(-pi/2,0,0);
        Rot2MF  = R_XZ2YZ*Rot2MF; 
end

% Generation of rotational matrix from molecular frame to lab frame
Rot2LF = R1_ZYZ_0(LF_Phi,LF_Psi,LF_Theta);

% Total rotation matrix
Rot_Total_XYZ = Rot2LF*Rot2MF;

% Coordination tranfrom of XYZ from mol. frame to lab frame
XYZ_T_R = (Rot_Total_XYZ*XYZ_T')';

%% Rotate the TDV and Raman tensor accordingly
% Determine if doing bond rotational average around z axis of mol. frame
if BondAvg
   R_Bond_Avg   = R1_ZYZ_1(0,0);
   R_Bond_Avg_2 = R2_ZYZ_1(0,0);
else
   R_Bond_Avg   = eye(3); 
   R_Bond_Avg_2 = eye(9);
end

% Total rotation matrix, include bond average
Rot_Total_mode = Rot2LF*R_Bond_Avg*Rot2MF;

% Coordination tranfrom of TDV from mol. frame to lab frame
TDV_Rot = (Rot_Total_mode*TDV')';

% Coordination tranfrom of Raman tensor from mol. frame to lab frame
% Raman_Rot = zeros(size(Raman));
% for i = 1:Mode_Num
%     Raman_Rot(:,:,i) = Rot_Total_mode*squeeze(Raman(:,:,i))*Rot_Total_mode';
% end

% kron(Aij,Bkl) = Cik,jl
Rot2MF_2 = kron(Rot2MF,Rot2MF);
Rot2LF_2 = kron(Rot2LF,Rot2LF);
Rot_Total_mode_2 = Rot2LF_2*R_Bond_Avg_2*Rot2MF_2;

Raman_V = reshape(Raman,9,[]);
Raman_Rot_V = Rot_Total_mode_2*Raman_V;
Raman_Rot_V(abs(Raman_Rot_V)<1E-10) = 0; % remove numerical error

Raman_Rot_M = reshape(Raman_Rot_V,3,3,[]);

%% Output
Output.LF            = S_Info;
Output.LF.Center_Ind = Center_Ind;
Output.LF.XYZ        = XYZ_T_R;
Output.LF.TDV        = TDV_Rot;
Output.LF.Raman_M    = Raman_Rot_M;
Output.LF.Raman_V    = Raman_Rot_V;

Output.Orig = S_Info;

Output.Transformation.C_G09_Frame  = C_G09_Frame;
Output.Transformation.Rot2MF       = Rot2MF;
Output.Transformation.Rot2LF       = Rot2LF;
Output.Transformation.Rot2MF_2     = Rot2MF_2;
Output.Transformation.Rot2LF_2     = Rot2LF_2;
Output.Transformation.R_Bond_Avg   = R_Bond_Avg;
Output.Transformation.R_Bond_Avg_2 = R_Bond_Avg_2;