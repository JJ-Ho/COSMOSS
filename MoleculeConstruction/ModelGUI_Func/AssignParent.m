function AssignParent(app,ParentInfo)
switch ParentInfo{1}
    case 'COSMOSS' %% Need to clean up UUID field later!
        app.ParentGUI    = ParentInfo{2};
        app.UIFigure.Tag = ParentInfo{2}.UUID;
        app.fileChooser  = ParentInfo{2}.fileChooser; % use the selected fileChooser in COSMOSS
        disp('Running sub-GUI from COSMOSS...')
    case 'Comb2'
        app.ParentGUI    = ParentInfo{2};
        app.UIFigure.Tag = ParentInfo{2}.UUID;
        app.fileChooser  = ParentInfo{2}.fileChooser; % use the selected fileChooser in COSMOSS
        app.CombOrder    = ParentInfo{3};
        disp('Running sub-GUI from Comb2...')
    case 'StandAlong'
        app.ParentGUI         = struct;
        app.ParentGUI.GUI_Tag = 'StandAlong';
        app.fileChooser = DefaultFileChooser;
        disp('Running Model GUI in stand alone mode...')
    otherwise
        disp([ParentInfo{1}, ' is not defined...'])
end