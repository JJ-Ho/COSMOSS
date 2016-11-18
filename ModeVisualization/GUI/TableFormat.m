function CF=TableFormat(SpecType)
% This function determine table format depends on the Spectrum type input

%% Defining column formats using ColumnFormat class
Index  = ColumnFormat('Index'  ,'short',40);
Freq   = ColumnFormat('Freq.'  ,'bank' ,50);
Norm1D = ColumnFormat('Norm 1D','bank' ,60);
Norm2D = ColumnFormat('Norm 2D','bank' ,60);
Mu_Int = ColumnFormat('Mu Int' ,'bank' ,50);
Mu_x   = ColumnFormat('Mu x'   ,'short',60);
Mu_y   = ColumnFormat('Mu y'   ,'short',60);
Mu_z   = ColumnFormat('Mu z'   ,'short',60);
A_Int  = ColumnFormat('Norm[A]','bank' ,50);
A_p1   = ColumnFormat('A p1'   ,'short',60);
A_p2   = ColumnFormat('A p2'   ,'short',60);
A_p3   = ColumnFormat('A p3'   ,'short',60);

%%
switch SpecType
    case 'FTIR'
        CF = merge2cell(Index,Freq,...
                        Mu_Int,Mu_x,Mu_y,Mu_z);
    case 'SFG'
        CF = merge2cell(Index,Freq,...
                        Norm1D,Norm2D,...
                        Mu_Int,Mu_x,Mu_y,Mu_z,...
                        A_Int,A_p1,A_p2,A_p3);
end