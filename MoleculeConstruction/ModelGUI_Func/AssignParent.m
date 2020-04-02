function AssignParent(app,ParentInfo)
switch ParentInfo{1}
    case 'COSMOSS'
        app.Parent      = ParentInfo{2};
        app.fileChooser = ParentInfo{2}.fileChooser; % use the selected fileChooser in COSMOSS
        disp('Running sub-GUI from COSMOSS...')
    case 'Comb2'
        app.Parent    = ParentInfo{2};
        app.fileChooser = ParentInfo{2}.fileChooser; % use the selected fileChooser in COSMOSS
        app.CombOrder = ParentInfo{3};
        disp('Running sub-GUI from Comb2...')
    case 'StandAlong'
        app.Parent.GUI_Tag = 'None';
        app.fileChooser    = DefaultFileChooser;
        disp('Running Model GUI in stand alone mode...')
end