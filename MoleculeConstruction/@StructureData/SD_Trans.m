function obj_T = SD_Trans(obj,V)
obj_T = SD_Copy(obj);

XYZ       = obj.XYZ;
LocCenter = obj.LocCenter;

% Make sure V is a row
if iscolumn(V)
    V = V';
end

% check if V is [1x3]
S = size(V);
if ~isequal(S,[1,3])
    disp('The size of V should be [1,3]...')
    disp('StructureData is not modified...')
    return
else
    XYZ_T       = bsxfun(@plus,XYZ,V);
    LocCenter_T = bsxfun(@plus,LocCenter,V);
        
    obj_T.XYZ       = XYZ_T;
    obj_T.LocCenter = LocCenter_T;
end   
