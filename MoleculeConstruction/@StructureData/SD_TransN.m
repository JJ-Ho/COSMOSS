function obj_TN = SD_TransN(obj,V,N)
% Given a translational vector and Number of copies N, this method duplicate 
% the orign StructureData by N times and move them by increatment of V.
%% Check Input
% Make sure V is a column
if isrow(V)
    V = V';
end

% check if V is [3x1]
if ~isequal(size(V),[3,1])
    disp('The size of V should be [1,3]...')
    disp('StructureData is not modified...')
    return
end

%% Assign properties that is not affected by either translation and duplication 
obj_TN = StructureData;
obj_TN.FilesName  = obj.FilesName;
obj_TN.hPlotFunc  = obj.hPlotFunc;
obj_TN.GUI_Inputs = obj.GUI_Inputs;

%% Duplicate propeties that are not affected by translational movement
obj_TN.AtomName  = repmat(obj.AtomName ,N,1);
obj_TN.LocFreq   = repmat(obj.LocFreq  ,N,1);
obj_TN.LocAnharm = repmat(obj.LocAnharm,N,1);
obj_TN.LocMu     = repmat(obj.LocMu    ,N,1);
obj_TN.LocAlpha  = repmat(obj.LocAlpha ,N,1);

%% Deal with properties that will be affected by the translational movement
V_T = reshape(bsxfun(@times,V,0:N-1), 1,3,[]);

XYZ_TN_3D  = bsxfun(@plus,repmat(obj.XYZ      ,1,1,N),V_T); 
obj_TN.XYZ = reshape(permute(      XYZ_TN_3D,[1,3,2]),[],3);

if ~isempty(obj.LocCenter)
    LocCenter_TN_3D  = bsxfun(@plus,repmat(obj.LocCenter,1,1,N),V_T);
    obj_TN.LocCenter = reshape(permute(LocCenter_TN_3D,[1,3,2]),[],3);
end