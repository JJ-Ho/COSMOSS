function S_PDB = Load_PDB(app,GUI_Inputs)
% This function read a PDB file and turn it into a StructureData object
% Copyright Jia-Jung Ho, 2013-2020

%% Get pdb file location 
[FilesName,PathName,~] = app.fileChooser.chooseFile('*.*');
                                
%% Parse molecule structure
if GUI_Inputs.Preprocessed
    % read the pre-processed MD sanpshots
    TextPattern = 'ATOM %*f %s %*s %*s %*f %f %f %f %*f %*f %*s';
    if iscell(FilesName)
        N_File = size(FilesName,2);
        fid = fopen([PathName FilesName{1}]);
        ATest = textscan(fid,TextPattern,'CollectOutput',1);
        fclose(fid);
        XYZ = zeros([size(ATest{2}),N_File]);
        
        progressbar;
        for i=1:N_File
            progressbar(i/N_File)
            fid = fopen([PathName FilesName{i}]);
            A = textscan(fid,TextPattern,'CollectOutput',1);
            fclose(fid);
            XYZ(:,:,i) = A{2};
        end
        
        AtomName  = A{1};
        C = strsplit(PathName,'/');
        FilesName  = C{end-1};
        
        % update GUI inputs in COSMOSS
        if ~isempty(app.Parent)
            app.Parent.CheckBox_Sampling.Value = 1;
            app.Parent.EditField_SampleN.Value = N_File;
            app.Parent.EditField_DD.Value = 0;
            app.Parent.EditField_ODD.Value = 0;
        end
        
    else
        fid = fopen([PathName FilesName]);
        A = textscan(fid,TextPattern,'CollectOutput',1);
        fclose(fid);
        
        AtomName   = A{1};
        XYZ(:,:,1) = A{2};
    end
else 
    N_Model = 1; % if there are more than one, select the first model in the PDB file by default
    PDB_orig = pdbread([PathName FilesName]);
    Atom_Data = PDB_orig.Model(N_Model).Atom;
    
    AtomName = {Atom_Data.AtomName}';
    XYZ      = [[Atom_Data.X]',[Atom_Data.Y]',[Atom_Data.Z]'];
    ChainID  = {Atom_Data.chainID};
end

%% output
S_PDB = StructureData;
S_PDB.XYZ       = XYZ;
S_PDB.AtomName  = AtomName;
S_PDB.FilesName = FilesName;
S_PDB.Extra.PDB_orig  = PDB_orig;
S_PDB.Extra.ChainID   = ChainID;

