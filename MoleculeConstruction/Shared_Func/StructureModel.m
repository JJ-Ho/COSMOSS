function [hStructure, ModelList, hPlotFunc] = StructureModel(StructModel,Handle)

%% List 
ModelList = {'1:Two Coupled Oscillators',...
             '2:PDB_AmideI',...
             '3:2D Grid',...
             '4:Ideal Betasheet',...
             '5:Combination of any two',...
             };
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
        hPlotFunc  = 'Non';
    case 1
        hStructure = Model_TCO(Export_handle);
        hPlotFunc  = 'PlotXYZfiles_Acid';
    case 2 
        hStructure = Model_PDB_AmideI(Export_handle);
        hPlotFunc  = 'PlotXYZfiles_AmideI';
    case 3
        hStructure = Model_TwoDGrid(Export_handle);
        hPlotFunc  = 'PlotXYZ_Grid';
    case 4
        hStructure = Model_Betasheet_AmideI(Export_handle);
        hPlotFunc  = 'Plot_Betasheet_AmideI';
    case 5
        hStructure = Model_Comb2(Export_handle);
        hPlotFunc  = 'Non';
    otherwise
        hStructure = 'Non';
        hPlotFunc  = 'Non';
        
        disp('Model List')
        disp('--------------------------')
        disp(ModelList')
        disp('Please select the models above...')
end
