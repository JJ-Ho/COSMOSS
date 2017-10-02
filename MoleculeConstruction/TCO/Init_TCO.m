function Init_TCO(app)
%% Initial vlues of the Orientation inputs
Orientation = zeros(3,7);
Orientation(1:2,1) = true;
Orientation(2,5) = 5;

ColumnName     = {       '', 'Phi', 'Psi','Theta',   'X',   'Y',   'Z'};
ColumnFormat   = {'logical','bank','bank', 'bank','bank','bank','bank'};
ColumnEditable = [     true,  true,  true,   true,  true,  true,  true];
ColumnWidth    = {       20,    30,    30,    40,    30,    30,    30};

app.UITable.Data           = Orientation;
app.UITable.ColumnName     = ColumnName;
app.UITable.ColumnFormat   = ColumnFormat;
app.UITable.ColumnEditable = ColumnEditable;
app.UITable.ColumnWidth    = ColumnWidth;
app.UITable.RowName        = 'numbered';