function obj_T = SD_Trans(obj,V)
    XYZ       = obj.XYZ;
    LocCenter = obj.LocCenter;

    % check if V is [1x3]
    S = size(V);
    if ~isequal(S,[1,3])
        if isequal(S,[3,1])
            disp('Transpost the input V from [3x1] to [1x3]...')
            V = V';
        else
            disp('The size of V should be [1,3]...')
            disp('StructureData is not modified...')
            return
        end
    else
        XYZ_T       = bsxfun(@plus,XYZ,V);
        LocCenter_T = bsxfun(@plus,LocCenter,V);

        obj_T           = obj;
        obj_T.XYZ       = XYZ_T;
        obj_T.LocCenter = LocCenter_T;
    end   
end