function obj_comb2 = SD_Comb2(obj_SD1,obj_SD2)
% Given two StructureDatas, this method will combie the two as a new
% StructureData
%% General properties
XYZ_1             = obj_SD1.XYZ;
AtomName_1        = obj_SD1.AtomName;
LocCenter_1       = obj_SD1.LocCenter;
LocFreq_1         = obj_SD1.LocFreq;
LocAnharm_1       = obj_SD1.LocAnharm;
DiagDisorder_1    = obj_SD1.DiagDisorder;
OffDiagDisorder_1 = obj_SD1.OffDiagDisorder;
LocMu_1           = obj_SD1.LocMu;
LocAlpha_1        = obj_SD1.LocAlpha;
Scaled_LocMu_1    = obj_SD1.Scaled_LocMu;
Scaled_LocAlpha_1 = obj_SD1.Scaled_LocAlpha;

XYZ_2             = obj_SD2.XYZ;
AtomName_2        = obj_SD2.AtomName;
LocCenter_2       = obj_SD2.LocCenter;
LocFreq_2         = obj_SD2.LocFreq;
LocAnharm_2       = obj_SD2.LocAnharm;
DiagDisorder_2    = obj_SD2.DiagDisorder;
OffDiagDisorder_2 = obj_SD2.OffDiagDisorder;
LocMu_2           = obj_SD2.LocMu;
LocAlpha_2        = obj_SD2.LocAlpha;
Scaled_LocMu_2    = obj_SD2.Scaled_LocMu;
Scaled_LocAlpha_2 = obj_SD2.Scaled_LocAlpha;

XYZ             = [            XYZ_1;           XYZ_2];
AtomName        = [       AtomName_1;       AtomName_2];
LocCenter       = [      LocCenter_1;      LocCenter_2];
LocFreq         = [        LocFreq_1;        LocFreq_2];
LocAnharm       = [      LocAnharm_1;      LocAnharm_2];
DiagDisorder    = [   DiagDisorder_1;   DiagDisorder_2];
OffDiagDisorder = [OffDiagDisorder_1;OffDiagDisorder_2];
LocMu           = [          LocMu_1;          LocMu_2];
LocAlpha        = [       LocAlpha_1;       LocAlpha_2];
Scaled_LocMu    = [   Scaled_LocMu_1;   Scaled_LocMu_2];
Scaled_LocAlpha = [Scaled_LocAlpha_1;Scaled_LocAlpha_2];

FilesName_1 = obj_SD1.FilesName;
FilesName_2 = obj_SD2.FilesName;
FilesName   = ['Comb2: ',FilesName_1,' & ',FilesName_2];

%% Deal with general extra properties
Extra_1 = obj_SD1.Extra;
Extra_2 = obj_SD2.Extra;
Extra.Extra_1 = Extra_1;
Extra.Extra_2 = Extra_2;

%% Deal with extra properties:AmideIAtomSerNo for peptides
N_Mode_1 = obj_SD1.Nmodes;
N_Mode_2 = obj_SD2.Nmodes;
N_Mode_total = N_Mode_1 + N_Mode_2;

AmideModeIndex = zeros(size(LocFreq,1),1);

if isfield(Extra_1,'AmideIAtomSerNo')
    AmideIAtomSerNo_1 = Extra_1.AmideIAtomSerNo;
    AmideModeIndex(1:N_Mode_1) = 1;
else
    AmideIAtomSerNo_1 =[];
end 

if isfield(Extra_2,'AmideIAtomSerNo')
    AmideIAtomSerNo_2 = Extra_2.AmideIAtomSerNo + obj_SD1.NAtoms;
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
Children = [obj_SD1,obj_SD2];

%% Reassign StructureData
obj_comb2 = StructureData;
obj_comb2.XYZ             = XYZ;
obj_comb2.AtomName        = AtomName;
obj_comb2.LocCenter       = LocCenter;
obj_comb2.LocFreq         = LocFreq;
obj_comb2.LocAnharm       = LocAnharm;
obj_comb2.DiagDisorder    = DiagDisorder;
obj_comb2.OffDiagDisorder = OffDiagDisorder;
obj_comb2.LocMu           = LocMu;
obj_comb2.LocAlpha        = LocAlpha;
obj_comb2.FilesName       = FilesName;
obj_comb2.Extra           = Extra;
obj_comb2.Children        = Children;
obj_comb2.Scaled_LocMu    = Scaled_LocMu;    % reassign to replace the automatically generated on from "obj_comb2.LocMu"
obj_comb2.Scaled_LocAlpha = Scaled_LocAlpha; % reassign to replace the automatically generated on from "obj_comb2.LocAlpha"
    