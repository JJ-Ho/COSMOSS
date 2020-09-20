function Export2Parent(app)
% This function push Structure data to COSMOSS after the Model strueture data updated

%% Decide if push structure to parent
switch app.ParentGUI.GUI_Tag
    case 'COSMOSS'
        app.ParentGUI.Structure = app.Structure;
        app.ParentGUI.UIFigure.Name = ['COSMOSS: ', app.UIFigure.Name]; % change Name of Main GUI to identify which Structural Model is using
        app.ParentGUI.RefreshOn % reset the all refresh tag to 1
    case 'Comb2'
        app.ParentGUI.RefreshSD % auto refresh the combined structure
    case 'StandAlong'
        % do nothing
    otherwise 
        disp('The ParentGUI type is not supported...')
end

%% Check which type of structure model is used and modify the main GUI accordingly
switch app.GUI_Tag
    case 'BSheet'
        app.ParentGUI.DropDown_ViewSamplingMode.Enable = 1;
        app.ParentGUI.isotopeDilutionCheckBox.Enable = 1;
        app.ParentGUI.isotopeDilutionCheckBox.Value  = 0;
    otherwise
        app.ParentGUI.DropDown_ViewSamplingMode.Enable = 0;
        app.ParentGUI.DropDown_ViewSamplingMode.Value  = app.ParentGUI.DropDown_ViewSamplingMode.Items(1);
        app.ParentGUI.isotopeDilutionCheckBox.Enable = 0;
        app.ParentGUI.isotopeDilutionCheckBox.Value  = 0;
end
    

%% Display update message
timeStamp    = datetime('now','TimeZone','local');
timeSamepStr = datestr(timeStamp,'yy/mm/dd HH:MM:SS');
disp([timeSamepStr,': ',app.GUI_Tag, ' Structure refreshed...'])
