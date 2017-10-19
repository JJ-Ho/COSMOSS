function ExportGUIData(app)

Export.Parent     = app.Parent;
Export.Structure  = app.Structure;
Export.CombOrder  = app.CombOrder;
Export.GUI_Inputs = app.Parse_GUI;

ExportName = ['Data_', app.GUI_Tag];
assignin('base',ExportName,Export)
disp(['Updated',ExportName,' exported!'])