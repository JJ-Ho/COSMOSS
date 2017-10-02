function GUI_Setup(app)
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
app.UITable_1D.ColumnWidth = {55,55,55};
app.UITable_1D.RowName = {'Ang.';'Pol.'};

Exp_2D = zeros(2,5);
Exp_2D(1,:) = 90;
app.UITable_2D.Data = Exp_2D;
app.UITable_2D.ColumnWidth = {44,44,40,41,40};
app.UITable_2D.RowName = {'Ang.';'Pol.'};

% Setup Lineshape list
[~,AvalibleL] = Conv_LineShape('List','','','');
app.DropDown_LineShape.Items = AvalibleL;

% Setup color map list
[~,CMap_List] = SelectColormap('List');
app.DropDown_ColorMap.Items = CMap_List;
