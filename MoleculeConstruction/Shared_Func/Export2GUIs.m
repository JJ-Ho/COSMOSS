function Export2GUIs(app)
% This function push Structure data to COSMOSS after the Model strueture data updated
%% check if the upper layer is COSMOSS, if yes push Structure to COSMOSS
switch app.Parent.GUI_Tag
    case 'COSMOSS'
        app.Parent.Structure = app.Structure;
        app.Parent.UIFigure.Name = ['COSMOSS: ', app.UIFigure.Name]; % change Name of Main GUI to identify which Structural Model is using
        app.Parent.RefreshOn % reset the all refresh tag to 1
    case 'Comb2'
        
end

