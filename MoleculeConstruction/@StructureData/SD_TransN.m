function obj_TN = SD_TransN(obj,V,N)
obj_TN = StructureData;

%% Assign properties that is not affected by either translation and duplication 
obj_TN.FilesName     = obj.FilesName;
obj_TN.hPlotFunc     = obj.hPlotFunc;
obj_TN.hParseGUIFunc = obj.hParseGUIFunc;
obj_TN.hGUIs         = obj.hGUIs; 

% obj_TN.Children = obj; 
% obj_TN.StructModel = 3; 
% obj_TN.Extra = '';

%% Duplicate propeties that are not affected by translational movement
obj_TN.AtomName  = repmat(obj.AtomName ,N,1);
obj_TN.LocFreq   = repmat(obj.LocFreq  ,N,1);
obj_TN.LocAnharm = repmat(obj.LocAnharm,N,1);
obj_TN.LocMu     = repmat(obj.LocMu    ,N,1);
obj_TN.LocAlpha  = repmat(obj.LocAlpha ,N,1);

%% Deal with properties that will be affected by the translational movement
% turn translational vector into column
if isrow(V)
    V = V';
end
V_T = reshape(bsxfun(@times,V,0:N-1), 1,3,[]);

XYZ_TN_3D       = bsxfun(@plus,repmat(obj.XYZ      ,1,1,N),V_T); 
LocCenter_TN_3D = bsxfun(@plus,repmat(obj.LocCenter,1,1,N),V_T);

obj_TN.XYZ       = reshape(permute(      XYZ_TN_3D,[1,3,2]),[],3);
obj_TN.LocCenter = reshape(permute(LocCenter_TN_3D,[1,3,2]),[],3);