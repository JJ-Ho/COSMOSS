function obj_T = SD_Trans(obj_SD,V)
% Given a translational vector V, this method move the StructureData accordingly 
%% Check Input
% Make sure V is a row
if iscolumn(V)
    V = V';
end

% check if V is [1x3]
if ~isequal(size(V),[1,3])
    disp('The size of V should be [1,3]...')
    disp('StructureData is not modified...')
    return
end

%% Copy the original obj to a new obj
obj_T = SD_Copy(obj_SD);

%% Apply translation to relevent properties
obj_T.XYZ = bsxfun(@plus,obj_SD.XYZ,V);

% check if the center location of modes is assigned
if ~isempty(obj_SD.LocCenter)
    obj_T.LocCenter = bsxfun(@plus,obj_SD.LocCenter,V);
end
