function [hModel, ModelList, hPlotFunc] = StructureModel(StructModel)
%% List 
ModelList = {'1:Two Coupled Oscillators',...
             '2:Ideal Betasheet',...
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
    case 1
        %hModel     = Model_TCO(Export_handle);
        hModel     = @Model_TCO;
        hPlotFunc  = @PlotXYZfiles_Acid;
    case 2
        hModel     = @Model_Betasheet_AmideI;
        hPlotFunc  = @Plot_Betasheet_AmideI;
    otherwise
        hModel     = 'Non';
        hPlotFunc  = 'Non';
        
        disp('Model List')
        disp('--------------------------')
        disp(ModelList')
        disp('Please select the models above...')
end
