function Export2GUIs(GUI_data)
% This function push Structure data to COSMOSS or Comb2 after the Model 
% strueture data updated

%% check if this program run stand along, if not push Structure info to hMain
if isfield(GUI_data,'hCOSMOSS')
    GUI_data_COSMOSS           = guidata(GUI_data.hCOSMOSS);
    GUI_data_COSMOSS.Structure = GUI_data.Structure;
    guidata(GUI_data.hCOSMOSS,GUI_data_COSMOSS)
    
    % change Name of Main GUI to help identifying which Structural Model is
    % using
    Model_Name    = GUI_data.hModel.Name;
    GUI_data.hCOSMOSS.Name = ['COSMOSS: ' Model_Name];
end

%% if belong to comb2, push GUI handle to comb2
if isfield(GUI_data,'Comb2_Order')
    GUI_data_Comb2 = guidata(GUI_data.hModel_Comb2);
    
    switch GUI_data.Comb2_Order
        case 1
            GUI_data_Comb2.hStruc1        = GUI_data.hModel;
        case 2
            GUI_data_Comb2.hStruc2        = GUI_data.hModel;
    end
    
    guidata(GUI_data.hModel_Comb2,GUI_data_Comb2);
    
    % push comb2 order # to GUI title, if necessary
    TitleStr = GUI_data.hModel.Name;
    if ~strcmp(TitleStr(1),'#')
        GUI_data.hModel.Name = ['#',int2str(GUI_data.Comb2_Order),': ',TitleStr];
    end
end
