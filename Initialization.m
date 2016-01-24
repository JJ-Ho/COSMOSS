% Initialization Script

Current_path = pwd;

% Add path of SubFunctions
Functions = genpath([Current_path '/SpecFunc']);
addpath(Functions);

GUIs = genpath([Current_path '/GUI']);
addpath(GUIs);

MoleculeConstruction = genpath([Current_path '/MoleculeConstruction']);
addpath(MoleculeConstruction);

ModeVisualization = genpath([Current_path '/ModeVisualization']);
addpath(ModeVisualization);

addpath([Current_path '/StructureFiles']);