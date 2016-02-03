% Initialization Script
Current_path = pwd;
addpath(Current_path);

% Add path to SubFunctions
SpectalFunctions = genpath([Current_path '/SpecFunc']);
addpath(SpectalFunctions);

MoleculeConstruction = genpath([Current_path '/MoleculeConstruction']);
addpath(MoleculeConstruction);

addpath([Current_path '/StructureFiles']);

ModeVisualization = genpath([Current_path '/ModeVisualization']);
addpath(ModeVisualization);

Main_GUIs = genpath([Current_path '/Main_GUI']);
addpath(Main_GUIs);
