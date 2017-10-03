function Structure = Build_PDB_AmideI_Structure(PDB,GUI_Inputs)    
%% Construct molecule
Num_Atoms = PDB.Num_Atoms;
XYZ       = PDB.XYZ;
AtomName  = PDB.AtomName;
FilesName = PDB.FilesName;
N_File    = PDB.N_File;

% test # modes and pre-allocate matix
Tmp1 = GetAmideI(XYZ(:,:,1),AtomName,FilesName,GUI_Inputs);
Nmodes = Tmp1.Nmodes;
Tmp_LocMu     = zeros(Nmodes,3,N_File);
Tmp_LocAlpha  = zeros(Nmodes,9,N_File);
Tmp_LocCenter = zeros(Nmodes,3,N_File);
Tmp_XYZ       = zeros(Num_Atoms,3,N_File);

for i = 1:N_File
    Tmp = GetAmideI(XYZ(:,:,i),AtomName,FilesName,GUI_Inputs);
    Tmp_LocMu(:,:,i)     = Tmp.LocMu;
    Tmp_LocAlpha(:,:,i)  = Tmp.LocAlpha;
    Tmp_LocCenter(:,:,i) = Tmp.LocCenter;
    Tmp_XYZ(:,:,i)       = Tmp.XYZ;
end

Structure = StructureData;
Structure.XYZ       = Tmp_XYZ;
Structure.AtomName  = Tmp1.AtomName;

Structure.LocCenter = Tmp_LocCenter;
Structure.LocFreq   = Tmp1.LocFreq;
Structure.LocAnharm = Tmp1.LocAnharm;
Structure.LocMu     = Tmp_LocMu;
Structure.LocAlpha  = Tmp_LocAlpha;

Structure.FilesName = Tmp1.FilesName;
Structure.Extra.AmideIAtomSerNo = Tmp1.Extra.AmideIAtomSerNo;