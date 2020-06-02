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
end

%% Display update message
timeStamp    = datetime('now','TimeZone','local');
timeSamepStr = datestr(timeStamp,'yy/mm/dd HH:MM:SS');
disp([timeSamepStr,': ',app.GUI_Tag, ' Structure refreshed...'])
