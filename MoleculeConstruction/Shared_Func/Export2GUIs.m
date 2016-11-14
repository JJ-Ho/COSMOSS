function Export2GUIs(GUI_data)
% This function push Structure data to COSMOSS after the Model 
% strueture data updated

%% check if the upper layer is COSMOSS, if yes push Structure to COSMOSS
if isfield(GUI_data,'hCOSMOSS')
    GUI_data_COSMOSS           = guidata(GUI_data.hCOSMOSS);
    GUI_data_COSMOSS.Structure = GUI_data.Structure;
    guidata(GUI_data.hCOSMOSS,GUI_data_COSMOSS)
    
    % change Name of Main GUI to help identifying which Structural Model is
    % using
    Model_Name    = GUI_data.hModel.Name;
    GUI_data.hCOSMOSS.Name = ['COSMOSS: ' Model_Name];
end

