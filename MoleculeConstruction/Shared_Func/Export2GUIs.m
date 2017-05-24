function Export2GUIs(GUI_data)
% This function push Structure data to COSMOSS after the Model 
% strueture data updated

%% check if the upper layer is COSMOSS, if yes push Structure to COSMOSS
if isfield(GUI_data,'hCOSMOSS')
    GUI_data_COSMOSS           = guidata(GUI_data.hCOSMOSS);
    GUI_data_COSMOSS.Structure = GUI_data.Structure;
    
    % reset the all refresh tag to 1
    Fname = fieldnames(GUI_data_COSMOSS.Refresh);
    for i = 1:length(Fname)
        GUI_data_COSMOSS.Refresh.(Fname{i}) = 1;
    end
    
    guidata(GUI_data.hCOSMOSS,GUI_data_COSMOSS)
    
    % change Name of Main GUI to help identifying which Structural Model is
    % using
    Model_Name    = GUI_data.hModel.Name;
    GUI_data.hCOSMOSS.Name = ['COSMOSS: ' Model_Name];
end

