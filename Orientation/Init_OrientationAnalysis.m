function Init_OrientationAnalysis(app)
%% Initial vlues of the Orientation inputs
Data = '';

ColumnName     = {       '', 'Intensity', 'Ratio'};
ColumnFormat   = {'logical',      'char',  'bank'};
ColumnEditable = [     true,       false,    true];
ColumnWidth    = {       20,         80,      60};

app.UITable.Data           = Data;
app.UITable.ColumnName     = ColumnName;
app.UITable.ColumnFormat   = ColumnFormat;
app.UITable.ColumnEditable = ColumnEditable;
app.UITable.ColumnWidth    = ColumnWidth;
app.UITable.RowName        = 'numbered';