function BuildModel(GUI_data, Comb2_Order)
%% Construct model GUI 
if isfield(GUI_data,'Load')
    StructModel = GUI_data.Load.StructModel;
else
    hGUIs       = GUI_data.hGUIs;
    switch Comb2_Order
        case 1
            StructModel = get(hGUIs.StructBox1,'Value');
        case 2
            StructModel = get(hGUIs.StructBox2,'Value');
    end
end 

[fhStructure,~,~] = StructureModel(StructModel);

hModel = feval(fhStructure,'Comb2',GUI_data.hModel_Comb2,Comb2_Order);

% Pass hCOSMOSS to the fig file of each sub-GUIs so they can always access
% COSMOSS if needed
hModel.UserData = GUI_data.hModel_Comb2.UserData;

%% Update GUI_data_Model
% pass empty Structure for ???
GUI_data_Model              = guidata(hModel);
GUI_data_Model.Structure    = [];

% Update the GUI inputs if Load Structure 
if isfield(GUI_data,'Load')
    hGUIs_Model = GUI_data_Model.hGUIs;
    Para_S      = GUI_data.Load.GUI_Inputs;
    FieldName   = GUI_data.Load.GUI_FieldName;
    
    UpdateGUIs(hGUIs_Model,Para_S,FieldName)
    
    % push the loaded structure data to sub-GUI
    GUI_data_Model.GUI_FieldName = FieldName;
    GUI_data_Model.Structure     = GUI_data.Load.Structure;
end

guidata(hModel,GUI_data_Model)

%% Update GUI_Data in Comb2 and propergate to COSMOSS
% Save the modle handle to proper name
switch Comb2_Order
    case 1
        GUI_data.hStruc1 = hModel;
    case 2
        GUI_data.hStruc2 = hModel;
end

guidata(GUI_data.hModel_Comb2,GUI_data)

Export2GUIs(GUI_data_Model);
