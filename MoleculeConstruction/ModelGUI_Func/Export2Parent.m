function Export2Parent(app)
% This function push Structure data to COSMOSS after the Model strueture data updated

%% Decide if push structure to parent
switch app.Parent.GUI_Tag
    case 'COSMOSS'
        app.Parent.Structure = app.Structure;
        app.Parent.UIFigure.Name = ['COSMOSS: ', app.UIFigure.Name]; % change Name of Main GUI to identify which Structural Model is using
        app.Parent.RefreshOn % reset the all refresh tag to 1
    case 'Comb2'
        app.Parent.RefreshSD % auto refresh the combined structure
    case 'StandAlong'
        % do nothing
end

%% Display update message
timeStamp    = datetime('now','TimeZone','local');
timeSamepStr = datestr(timeStamp,'yy/mm/dd HH:MM:SS');
disp([timeSamepStr,': ',app.GUI_Tag, ' Structure refreshed...'])
