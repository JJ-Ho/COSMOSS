
function [hModel, ModelList, hPlotFunc, hGUIParser] = StructureModel(StructModel)
%% List 
ModelList = {'1:Two coupled oscillators',...
             '2:Extract Amide-I from PDB files',...
             '3:Build 2D grid from G09 monomer',...
             '4:Ideal betasheet',...
             '5:Combination of any two above',...
             };
%% Run Models

% check if input is an figure handle
% if ishandle(Handle)
%     Export_handle = Handle;
% else
%     Export_handle = 'StandAlone';
% end

switch StructModel
    case 0 % for exporting ModelList only
        hModel     = 'Non'; 
        hPlotFunc  = 'Non';
        hGUIParser = 'Non';
    case 1
        hModel     = @Model_TCO;
        hPlotFunc  = @PlotXYZfiles_Acid;
        hGUIParser = @ParseGUI_TCO;
    case 2 
        hModel     = @Model_PDB_AmideI;
        hPlotFunc  = @PlotXYZfiles_AmideI;
        hGUIParser = @ParseGUI_AmideI;
    case 3
        hModel     = @Model_TwoDGrid;
        hPlotFunc  = @PlotXYZ_Grid;
        hGUIParser = @ParseGUI_TwoDGrid;
    case 4
        hModel     = @Model_Betasheet_AmideI;
        hPlotFunc  = @Plot_Betasheet_AmideI;
        hGUIParser = @ParseGUI_Betasheet;
    case 5
        hModel     = @Model_Comb2;
        hPlotFunc  = 'Non';
        hGUIParser = @ParseGUI_Comb2;
    otherwise
        hModel     = 'Non';
        hPlotFunc  = 'Non';
        hGUIParser = 'Non';
        
        disp('Model List')
        disp('--------------------------')
        disp(ModelList')
        disp('Please select the models above...')
end
