function Output = LoadStructureData(ParentType,hParent,varargin)
PWD = pwd;
Default_folder = [PWD, '/StructureFiles/'];

[FilesName,PathName] = uigetfile({'*.mat','StructureData class'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',Default_folder);

load([PathName FilesName]);
if exist('SD','var')
    Structure = SD; % rename to previously saved StrcutureData (SD)
    disp([FilesName,' loaded...'])
else
    disp([FilesName,' does not have a StructureData named SD '])
    return
end

%% Check ParentType
switch ParentType
    case 'COSMOSS'
        Output = 'No GUI handel for the LoadStructureData case';
        hParent.Structure = Structure; % directly push the StructureData to COSMOSS
    case 'Comb2'
        Output.Structure  = Structure; % save the StructureData in a fake GUI handle "Output"
    otherwise
        Output = Structure; % diretly export the StructureData
end