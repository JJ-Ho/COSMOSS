function S_PDB = Load_PDB(app,GUI_Inputs)
% This function read a PDB file and turn it into a StructureData object
% Copyright Jia-Jung Ho, 2013-2020

%% Get pdb file location 
PWD = pwd;
PDB_Path = [PWD, '/StructureFiles/PDB/'];

[FilesName,PathName,~] = uigetfile({'*.pdb','PDB file'; ...
                                    '*,*','All Files'},...
                                    'MultiSelect','on',...
                                    'Select inputs',PDB_Path);
                                
%% Parse molecule structure
if GUI_Inputs.Preprocessed
    % read the pre-processed MD sanpshots with the same molecule
    TextPattern = 'ATOM %*f %s %*s %*s %*f %f %f %f %*f %*f %*s';
    if iscell(FilesName)
        N_File = size(FilesName,2);
        fid = fopen([PathName FilesName{1}]);
        ATest = textscan(fid,TextPattern,'CollectOutput',1);
        fclose(fid);
        XYZ = zeros([size(ATest{2}),N_File]); % Creat all zeros RawData Matrix
        
        progressbar;
        for i=1:N_File
            progressbar(i/N_File)
            % read preprocessed pdb file
            fid = fopen([PathName FilesName{i}]);
            A = textscan(fid,TextPattern,'CollectOutput',1);
            fclose(fid);

            XYZ(:,:,i) = A{2};
        end
        AtomName  = A{1};
        Num_Atoms = size(XYZ,1);
        
        C = strsplit(PathName,'/');
        FilesName  = C{end-1};
        
        % deal with the GUI inputs in COSMOSS
        if ~isempty(app.Parent)
            app.Parent.CheckBox_Sampling.Value = 1;
            app.Parent.EditField_SampleN.Value = N_File;
            app.Parent.EditField_DD.Value = 0;
            app.Parent.EditField_ODD.Value = 0;
        end
        
    else
        N_File = 1;
        fid = fopen([PathName FilesName]);
        A = textscan(fid,TextPattern,'CollectOutput',1);
        fclose(fid);
        
        AtomName   = A{1};
        XYZ(:,:,1) = A{2};
        Num_Atoms = size(XYZ,1);
    end
else 
    N_File = 1;
    PDB_orig = pdbread([PathName FilesName]);
    Atom_Data = PDB_orig.Model.Atom;
    Num_Atoms = size(Atom_Data,2);

    % Get coordination data
    XYZ = zeros(Num_Atoms,3);
    AtomName = cell(Num_Atoms,1);
    for II = 1:Num_Atoms
        A = Atom_Data(II);
        XYZ(II,:) = [A.X, A.Y, A.Z];
        AtomName{II} = Atom_Data(II).AtomName;
    end
end

%% output
S_PDB = StructureData;
S_PDB.XYZ       = XYZ;
S_PDB.AtomName  = AtomName;
S_PDB.FilesName = FilesName;
S_PDB.Extra.PDB = PDB_orig;
