function obj = SD_Rot(obj,Phi,Psi,Theta)
% Translate the CoM to origin
CoM_0  = obj.CoM;
obj = SD_Trans(obj,-CoM_0);

% Apply rotation
XYZ       = obj.XYZ;
LocCenter = obj.LocCenter;
LocMu     = obj.LocMu;
LocAlpha  = obj.LocAlpha;

R1 = R1_ZYZ_0(Phi,Psi,Theta);
R2 = R2_ZYZ_0(Phi,Psi,Theta);

XYZ_R       = (R1*      XYZ')';
LocCenter_R = (R1*LocCenter')';
LocMu_R     = (R1*    LocMu')';
LocAlpha_R  = (R2* LocAlpha')';

obj.XYZ       = XYZ_R;
obj.LocCenter = LocCenter_R;
obj.LocMu     = LocMu_R;
obj.LocAlpha  = LocAlpha_R;

% Translate the CoM back to where it was
obj = SD_Trans(obj,CoM_0);
