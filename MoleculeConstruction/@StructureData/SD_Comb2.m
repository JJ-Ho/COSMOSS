function obj_comb2 = SD_Comb2(obj1,obj2)
%% General properties
XYZ_1       = obj1.XYZ;
AtomName_1  = obj1.AtomName;
LocCenter_1 = obj1.LocCenter;
LocFreq_1   = obj1.LocFreq;
LocAnharm_1 = obj1.LocAnharm;
LocMu_1     = obj1.LocMu;
LocAlpha_1  = obj1.LocAlpha;
Scaled_LocMu_1    = obj1.Scaled_LocMu;
Scaled_LocAlpha_1 = obj1.Scaled_LocAlpha;

XYZ_2       = obj2.XYZ;
AtomName_2  = obj2.AtomName;
LocCenter_2 = obj2.LocCenter;
LocFreq_2   = obj2.LocFreq;
LocAnharm_2 = obj2.LocAnharm;
LocMu_2     = obj2.LocMu;
LocAlpha_2  = obj2.LocAlpha;
Scaled_LocMu_2    = obj2.Scaled_LocMu;
Scaled_LocAlpha_2 = obj2.Scaled_LocAlpha;

XYZ        = [      XYZ_1;       XYZ_2];
AtomName   = [ AtomName_1;  AtomName_2];
LocCenter  = [LocCenter_1; LocCenter_2];
LocFreq    = [  LocFreq_1;   LocFreq_2];
LocAnharm  = [LocAnharm_1; LocAnharm_2];
LocMu      = [    LocMu_1;     LocMu_2];
LocAlpha   = [ LocAlpha_1;  LocAlpha_2];
Scaled_LocMu    = [   Scaled_LocMu_1;   Scaled_LocMu_2];
Scaled_LocAlpha = [Scaled_LocAlpha_1;Scaled_LocAlpha_2];

FilesName_1 = obj1.FilesName;
FilesName_2 = obj2.FilesName;
FilesName   = ['Comb2: ',FilesName_1,' & ',FilesName_2];

%% Deal with general extra properties
Extra_1 = obj1.Extra;
Extra_2 = obj2.Extra;
Extra.Extra_1 = Extra_1;
Extra.Extra_2 = Extra_2;

%% Deal with extra properties:AmideIAtomSerNo for peptides
N_Mode_1 = obj1.Nmodes;
N_Mode_2 = obj2.Nmodes;
N_Mode_total = N_Mode_1 + N_Mode_2;

AmideModeIndex = zeros(size(LocFreq,1),1);

if isfield(Extra_1,'AmideIAtomSerNo')
    AmideIAtomSerNo_1 = Extra_1.AmideIAtomSerNo;
    AmideModeIndex(1:N_Mode_1) = 1;
else
    AmideIAtomSerNo_1 =[];
end 

if isfield(Extra_2,'AmideIAtomSerNo')
    AmideIAtomSerNo_2 = Extra_2.AmideIAtomSerNo + obj1.NAtoms;
    AmideModeIndex(N_Mode_1+1:N_Mode_total) = 1;
else
    AmideIAtomSerNo_2 =[];
end

AmideIAtomSerNo = [AmideIAtomSerNo_1;AmideIAtomSerNo_2];
if ~isempty(AmideIAtomSerNo)
    Extra.AmideIAtomSerNo = AmideIAtomSerNo;
    Extra.AmideModeIndex  = logical(AmideModeIndex);
end

%% Inherent all Children
% N_Child1 = length(obj1.Children);
% N_Child2 = length(obj2.Children);
% 
% if N_Child1
%     Children = obj1.Children;
% else
%     Children = obj1;
% end
% 
% if N_Child2
%     Children = [Children,obj2.Children];
% else
%     Children = [Children,obj2];
% end
    
Children = [obj1,obj2];

%% Reassign StructureData
obj_comb2 = StructureData;
obj_comb2.StructModel = 5;
obj_comb2.XYZ         = XYZ;
obj_comb2.AtomName    = AtomName;
obj_comb2.LocCenter   = LocCenter;
obj_comb2.LocFreq     = LocFreq;
obj_comb2.LocAnharm   = LocAnharm;
obj_comb2.LocMu       = LocMu;
obj_comb2.LocAlpha    = LocAlpha;
obj_comb2.FilesName   = FilesName;
obj_comb2.Extra       = Extra;
obj_comb2.Children    = Children;
obj_comb2.Scaled_LocMu    = Scaled_LocMu;    % reassign to replace the automatically generated on from "obj_comb2.LocMu"
obj_comb2.Scaled_LocAlpha = Scaled_LocAlpha; % reassign to replace the automatically generated on from "obj_comb2.LocAlpha"
    