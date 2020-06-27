function ExportGUIData(app)
switch app.GUI_Tag
    case 'COSMOSS'
        Export.Structure  = app.Structure;
        Export.Data_FTIR  = app.Data_FTIR;
        Export.Data_SFG   = app.Data_SFG;
        Export.Data_2DIR  = app.Data_2DIR;
        Export.Data_2DSFG = app.Data_2DSFG;            
        %Export.hCOSMOSS   = app;
        %Export.hModel     = app.hModel;
    otherwise
        Export.Structure  = app.Structure;
        Export.CombOrder  = app.CombOrder;
        Export.GUI_Inputs = app.Parse_GUI;
        %Export.Parent     = app.Parent;
        %Export.app        = app;
end

ExportName = ['Data_', app.GUI_Tag];
assignin('base',ExportName,Export)
disp([ExportName,' is exported to workspace...'])