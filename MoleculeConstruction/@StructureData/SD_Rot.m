function obj_RT = SD_Rot(obj_SD,R1)
% This method rotate the StructureData object about it's CoM
%% Check Input
% check if R is [3x3]
if ~isequal(size(R1),[3,3])
    disp('The size of Rotational matrix "R" should be [3,3]...')
    disp('StructureData is not modified...')
    return
end

%% Translate the CoM to origin
obj_T = SD_Trans(obj_SD,-obj_SD.CoM);

%% Apply rotation to corresponding properties
obj_R     = SD_Copy(obj_T);
obj_R.XYZ = (R1*obj_T.XYZ')';

% Make sure if the properties exist
if ~isempty(obj_T.LocCenter)
    obj_R.LocCenter = (R1 * obj_T.LocCenter')';
end

if ~isempty(obj_T.LocMu)
    obj_R.LocMu = (R1 * obj_T.LocMu')';
end

if ~isempty(obj_T.LocAlpha)
    R2 = kron(R1,R1);
    obj_R.LocAlpha = (R2 * obj_T.LocAlpha')';
end

%% Translate the CoM back to where it was
obj_RT = SD_Trans(obj_R,obj_SD.CoM);
