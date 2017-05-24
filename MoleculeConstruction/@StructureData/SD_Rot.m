function obj_R = SD_Rot(obj,Phi,Psi,Theta)
% Translate the CoM to origin
CoM_0  = obj.CoM;
obj_T0 = SD_Trans(obj,-CoM_0);

% Apply rotation
XYZ       = obj_T0.XYZ;
LocCenter = obj_T0.LocCenter;
LocMu     = obj_T0.LocMu;
LocAlpha  = obj_T0.LocAlpha;

R1 = R1_ZYZ_0(Phi,Psi,Theta);
R2 = R2_ZYZ_0(Phi,Psi,Theta);

XYZ_R       = (R1*      XYZ')';
LocCenter_R = (R1*LocCenter')';
LocMu_R     = (R1*    LocMu')';
LocAlpha_R  = (R2* LocAlpha')';

obj_T0.XYZ       = XYZ_R;
obj_T0.LocCenter = LocCenter_R;
obj_T0.LocMu     = LocMu_R;
obj_T0.LocAlpha  = LocAlpha_R;

% Translate the CoM back to where it was
obj_R = SD_Trans(obj_T0,CoM_0);
