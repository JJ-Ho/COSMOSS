function obj_RT = SD_Rot(obj,R1)
% Translate the CoM to origin
obj_T = SD_Trans(obj,-obj.CoM);

% Apply rotation
XYZ       = obj_T.XYZ;
LocCenter = obj_T.LocCenter;
LocMu     = obj_T.LocMu;
LocAlpha  = obj_T.LocAlpha;

R2 = kron(R1,R1);
XYZ_R       = (R1*      XYZ')';
LocCenter_R = (R1*LocCenter')';
LocMu_R     = (R1*    LocMu')';
LocAlpha_R  = (R2* LocAlpha')';

obj_R = SD_Copy(obj_T);
obj_R.XYZ       = XYZ_R;
obj_R.LocCenter = LocCenter_R;
obj_R.LocMu     = LocMu_R;
obj_R.LocAlpha  = LocAlpha_R;

% Translate the CoM back to where it was
obj_RT = SD_Trans(obj_R,obj.CoM);
