function [hModel, ModelList] = StructureModel(StructModel)
% This function is the indexing reference of all the structure construction
% model. Only COSMOSS and Comb2 use it to call model generation GUIs

%% List 
ModelList = {...
             '1:Two coupled oscillators',...
             '2:Extract Amide-I from PDB files',...
             '3:Build 2D grid from G09 monomer',...
             '4:Ideal betasheet',...
             '5:Load previously constructed model',...
             '6:Construct cavity modes',...
             '7:Combination of any two above',...
             };
         
%% Run Models
switch StructModel
    case 1
        hModel = @Model_TCO;
    case 2 
        hModel = @Model_PDB_AmideOne;
    case 3
        hModel = @Model_TwoDGrid;
    case 4
        hModel = @Model_Betasheet_AmideI;
    case 5
        hModel = @LoadStructureData;
    case 6
        hModel = @Model_CavityModes;
    case 7
        hModel = @Model_Comb2;

    case 'List'
        % silently output Model list
        hModel     = 'Non';
        
    otherwise
        hModel     = 'Non';
        disp('Model List')
        disp('--------------------------')
        disp(ModelList')
        disp('Please select the models above...')
end
