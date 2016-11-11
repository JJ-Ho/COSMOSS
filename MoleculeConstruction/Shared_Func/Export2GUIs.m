function Export2GUIs(handles)

% check if this program run stand along, if not push Structure info to hMain
if isfield(handles,'hMain')
    Data_Main = guidata(handles.hMain);
    Data_Main.Structure = handles.Structure;
    guidata(handles.hMain,Data_Main)
    
    % change Name of Main GUI to help identifying which Structural Model is
    % using
    Model_Name    = handles.hModel.Name;
    handles.hMain.Name = ['COSMOSS: ' Model_Name];
end

% if belong to comb2, push GUI handle to comb2
if isfield(handles,'hComb2')
    Data_Comb2 = guidata(handles.hComb2);
    
    % update model specific plotting function
    StructModel = handles.Structure.StructModel;
    [~,~,hPlotFunc,~] = StructureModel(StructModel);
        
    switch handles.Comb2_Order
        case 1
            Data_Comb2.hStruc1    = handles.hModel;
            Data_Comb2.hPlotFunc1 = hPlotFunc;
            Data_Comb2.GUI_FieldName1 = handles.GUI_FieldName;
        case 2
            Data_Comb2.hStruc2    = handles.hModel;
            Data_Comb2.hPlotFunc2 = hPlotFunc;
            Data_Comb2.GUI_FieldName2 = handles.GUI_FieldName;
    end
    
    guidata(handles.hComb2,Data_Comb2);
    
    % push comb2 order # to GUI title, if necessary
    TitleStr = handles.hModel.Name;
    if ~strcmp(TitleStr(1),'#')
        handles.hModel.Name = ['#',int2str(handles.Comb2_Order),': ',TitleStr];
    end
end
