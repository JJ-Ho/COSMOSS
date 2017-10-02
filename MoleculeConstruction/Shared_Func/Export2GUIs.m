function Export2GUIs(app)
% This function push Structure data to COSMOSS after the Model 
% strueture data updated

%% check if the upper layer is COSMOSS, if yes push Structure to COSMOSS
if ~isempty(app.Parent)
    
    app.Parent.hModel    = app;
    app.Parent.Structure = app.Structure;
    app.Parent.UIFigure.Name = ['COSMOSS: ', app.UIFigure.Name];     % change Name of Main GUI to help identifying which Structural Model is using
%     % reset the all refresh tag to 1
%     Fname = fieldnames(GUI_data_COSMOSS.Refresh);
%     for i = 1:length(Fname)
%         GUI_data_COSMOSS.Refresh.(Fname{i}) = 1;
%     end
%     
%     guidata(app.hCOSMOSS,GUI_data_COSMOSS)
    
end

