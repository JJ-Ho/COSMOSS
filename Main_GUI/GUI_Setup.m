function GUI_Setup(app)
% Setup initial values of the experimental setups
Exp_1D = zeros(2,3);
Exp_1D(1,:) = 90;
app.UITable_1D.Data = Exp_1D;
app.UITable_1D.ColumnWidth = {55,55,55};
app.UITable_1D.RowName = {'Ang.';'Pol.'};

Exp_2D = zeros(2,5);
Exp_2D(1,:) = 90;
app.UITable_2D.Data = Exp_2D;
app.UITable_2D.ColumnWidth = {44,44,40,41,40};
app.UITable_2D.RowName = {'Ang.';'Pol.'};