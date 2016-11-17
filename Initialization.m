% Initialization Script
Current_path = pwd;
addpath(Current_path);

% Add path to SubFunctions
SpectalFunctions = genpath([Current_path '/SpecFunc']);
addpath(SpectalFunctions);

MoleculeConstruction = genpath([Current_path '/MoleculeConstruction']);
addpath(MoleculeConstruction);

StructureFiles = genpath([Current_path '/StructureFiles']);
addpath(StructureFiles);

ModeVisualization = genpath([Current_path '/ModeVisualization']);
addpath(ModeVisualization);

Main_GUIs = genpath([Current_path '/Main_GUI']);
addpath(Main_GUIs);

ServerVersion = genpath([Current_path '/ServerVersion']);
addpath(ServerVersion);

% Add path to SubFunctions
SpectalFunctions = genpath([Current_path '/GenerationCode']);
addpath(SpectalFunctions);