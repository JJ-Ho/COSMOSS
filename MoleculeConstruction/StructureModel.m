function [hStructure, ModelList] = StructureModel(StructModel,Handle)

%% List 
ModelList = {'1:Two Coupled Oscillators',...
             '2:PDB_AmideI',...
             '3:2D Grid',...
             '4:Ideal Betasheet'};
%% Run Models

% check if input is an figure handle
if ishandle(Handle)
    Export_handle = Handle;
else
    Export_handle = 'StandAlone';
end

switch StructModel
    case 0
        hStructure = 'Non'; % for exporting ModelList only
    case 1
        hStructure = Model_TCO(Export_handle);
    case 2 
        hStructure = Model_PDB_AmideI(Export_handle);
    case 3
        hStructure = Model_TwoDGrid(Export_handle);
    case 4
        hStructure = Model_Betasheet_AmideI(Export_handle);
    otherwise
        hStructure = 'Non';
        
        disp('Model List')
        disp('--------------------------')
        disp(ModelList')
        disp('Please select the models above...')
end
