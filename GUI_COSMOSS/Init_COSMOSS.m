function Init_COSMOSS(app)
%% Setup paths
Current_path = pwd;
addpath(Current_path);

% Add path 
MoleculeConstruction = genpath([Current_path '/MoleculeConstruction']);
addpath(MoleculeConstruction);

%             StructureFiles = genpath([Current_path '/StructureFiles']);
%             addpath(StructureFiles);
%             
SpectalFunctions = genpath([Current_path '/SpecFunc']);
addpath(SpectalFunctions);
%             
%             ServerVersion = genpath([Current_path '/ServerVersion']);
%             addpath(ServerVersion);
%             
%             AnalysisTools = genpath([Current_path '/AnalysisTools']);
%             addpath(AnalysisTools);
%             
%             ModeVisualization = genpath([Current_path '/ModeVisualization']);
%             addpath(ModeVisualization);


%% Place default values on GUI
% Initialize refresh tags
RefreshOn(app)

% Setup the Model structure list
[~,StructureModelList,~] = StructureModel(0);
app.ListBox_Model.Items = StructureModelList;
app.ListBox_Model.ItemsData = 1:length(StructureModelList);

% Setup Coupling model items
[~,CouplingList] = Coupling('','List','');
app.DropDown_Coupling.Items = CouplingList;

% Setup Ensemble average list
[~,~,R_List,M_List] = LabFrameAvg('List','List','');
app.DropDown_RotAvg.Items = R_List;
app.DropDown_MirrorPlane.Items = M_List;

% Setup initial values of the experimental setups
Exp_1D = zeros(2,3);
Exp_1D(1,:) = 90;
app.UITable_1D.Data = Exp_1D;
app.UITable_1D.ColumnWidth = {40,44,40};
app.UITable_1D.ColumnEditable = true;
app.UITable_1D.RowName = {'Ang.';'Pol.'};

Exp_2D = zeros(2,5);
Exp_2D(1,:) = 90;
app.UITable_2D.Data = Exp_2D;
app.UITable_2D.ColumnWidth = {44,44,40,41,40};
app.UITable_2D.ColumnEditable = true;
app.UITable_2D.RowName = {'Ang.';'Pol.'};

% Setup Lineshape list
[~,AvalibleL] = Conv_LineShape('List','','','');
app.DropDown_LineShape.Items = AvalibleL;

% Setup color map list
[~,CMap_List] = SelectColormap('List');
app.DropDown_ColorMap.Items = CMap_List;

%% Add Properties Change Listner
SA = [...
    app.DropDown_LocFreq;...
    app.DropDown_Coupling;...
    app.EditField_NN;...
    app.EditField_FeynmannCutoff;...
    app.DropDown_RotAvg;...
    app.DropDown_MirrorPlane;...
    app.CheckBox_Sampling;...
    app.EditField_SampleN;...
    app.CheckBox_DynamicFigUpdate;...
    app.EditField_FluctuationCorrelation;...
    app.EditField_DD;...
    app.EditField_ODD;...
    ];

for i = 1:length(SA)
    SA(i).ValueChangedFcn = @(~,~)RefreshOn(app);
end

app.UITable_1D.CellEditCallback = @(~,~)RefreshOn(app);
app.UITable_2D.CellEditCallback = @(~,~)RefreshOn(app);
