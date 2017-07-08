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

%% propagate the action to Children
NChild = length(obj.Children);
if NChild
    for i = 1:NChild
        obj_Child = obj.Children(i);
        CoM_Child = obj_Child.CoM;

        % Translate the CoM to origin
        obj_Child = SD_Trans(obj_Child,-CoM_Child);
        
        % Apply rotation
        XYZ_Child       = obj_Child.XYZ;
        LocCenter_Child = obj_Child.LocCenter;
        LocMu_Child     = obj_Child.LocMu;
        LocAlpha_Child  = obj_Child.LocAlpha;

        obj_Child.XYZ       = (R1*      XYZ_Child')';
        obj_Child.LocCenter = (R1*LocCenter_Child')';
        obj_Child.LocMu     = (R1*    LocMu_Child')';
        obj_Child.LocAlpha  = (R2* LocAlpha_Child')';

        % Translate the CoM back to where it was
        SD_Trans(obj_Child,CoM_Child);
       
    end
end
