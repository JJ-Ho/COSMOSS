function obj_T = SD_Trans(obj,V)
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
obj_T = SD_Copy(obj);

%% Apply translation to relevent properties
XYZ_T     = bsxfun(@plus,obj.XYZ,V);
obj_T.XYZ = XYZ_T;

% check if the center location of modes is assigned
if ~isempty(obj.LocCenter)
    LocCenter_T = bsxfun(@plus,obj.LocCenter,V);
    obj_T.LocCenter = LocCenter_T;
end
