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

%% Copy properties from input
obj_TN = SD_Copy(obj);

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

%% Generate Hamiltonian after producing translational copie
% note: its possible that the L_Index is outside the range of current mode
% number since SD_TransN method can be used during the strutural model
% construction, wher the struture model is not yet complete. As a result, I
% decide not to generate Hamiltonian information for this method.
% obj_TN = SD_1ExH(obj_TN);