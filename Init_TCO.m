function Init_TCO(app)
%% Initial vlues of the Orientation inputs
Orientation = zeros(3,4);
Orientation(1:2,1) = true;

ColumnName     = {       '', 'Phi', 'Psi','Theta'};
ColumnFormat   = {'logical','bank','bank', 'bank'};
ColumnEditable = [     true,  true,  true,  true];
ColumnWidth    = {       20,    40,    40,    40};

app.UITable.Data           = Orientation;
app.UITable.RowName        = 'numbered';
app.UITable.ColumnName     = ColumnName;
app.UITable.ColumnFormat   = ColumnFormat;
app.UITable.ColumnEditable = ColumnEditable;
app.UITable.ColumnWidth    = ColumnWidth;