function Init_COSMOSS_GUI(app)
%% Setup default GUI appearences
% Setup the list of structure models
[~,StructureModelList] = StructureModel('List');
app.ListBox_Model.Items  = StructureModelList;
app.ListBox_Model.ItemsData = 1:length(StructureModelList);

% Setup the list of ensemble averages
[~,~,R_List,M_List] = LabFrameAvg('List','List','');
app.DropDown_RotAvg.Items = R_List;
app.DropDown_MirrorPlane.Items = M_List;

% Setup the default values for the experimental setup
Exp_1D = zeros(2,3);
Exp_1D(1,:) = 90;
app.UITable_1D.Data = Exp_1D;

Exp_2D = zeros(2,5);
Exp_2D(1,:) = 90;
app.UITable_2D.Data = Exp_2D;

% Setup the list of Lineshapes
[~,AvalibleL] = Conv_LineShape('List','','','');
app.DropDown_LineShape.Items = AvalibleL;

% Setup the list of color maps
[~,CMap_List] = SelectColormap('List');
app.DropDown_ColorMap.Items = CMap_List;

%% Attach properties change listners
SA = [...
      app.EditField_FeynmannCutoff;...
      app.DropDown_RotAvg;...
      app.DropDown_MirrorPlane;...
      app.CheckBox_Sampling;...
      app.EditField_SampleN;...
      app.CheckBox_DynamicFigUpdate;...
      app.EditField_FluctuationCorrelation;...
      app.EditField_ODD;...
      app.UITable_1D;...
      app.UITable_2D;...
      ];

for i = 1:length(SA)
    if isprop(SA(i),'ValueChangedFcn')
        SA(i).ValueChangedFcn = @(~,~)RefreshOn(app);
        
    elseif isprop(SA(i),'CellEditCallback')
        SA(i).CellEditCallback = @(~,~)RefreshOn(app);
        
    else
        disp(['Cannot attach propertity listner to the object: ',SA(i)])
        
    end 
end

%% MISC tasks
app.GUI_Tag = 'COSMOSS'; % add tag for separating cases in ExportGUIData()