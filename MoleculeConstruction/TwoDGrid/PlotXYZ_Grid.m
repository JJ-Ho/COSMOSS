function hF = PlotXYZ_Grid(Structure)

%% Decide connectivity
Num_Modes = Structure.Num_Modes;
Atom_Num  = Structure.Monomer.Atom_Num;
XYZ       = Structure.XYZ;
Conn      = Connectivity(XYZ);

% Add S-C connection for Ester
for k = 1:Num_Modes
    Conn(13+(k-1)*Atom_Num,3+(k-1)*Atom_Num) = 1; 
end
%% Figures
hF = figure;
hold on
gplot3(Conn,XYZ);
hold off

%% Figure options
axis equal;
rotate3d on

xlabel('X')
ylabel('Y')
