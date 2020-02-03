function Update_OrientTable(app)
%% Reorginze Data_Brush
Data_Table = reshape([app.Data_Brush(:).TableRow],3,[])';
Row_Color = reshape([app.Data_Brush(:).Color],3,[])';

%% Initial vlues of the Orientation inputs
ColumnName     = {       '', 'Intensity', 'Ratio'};
ColumnFormat   = {'logical',      'char',  'bank'};
ColumnEditable = [     true,       false,    true];
ColumnWidth    = {       20,         80,      60};

app.UITable.ColumnName      = ColumnName;
app.UITable.ColumnFormat    = ColumnFormat;
app.UITable.ColumnEditable  = ColumnEditable;
app.UITable.ColumnWidth     = ColumnWidth;
app.UITable.RowName         = 'numbered';
app.UITable.BackgroundColor = Row_Color;
app.UITable.Data            = Data_Table;
