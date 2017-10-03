function O = Parse_TCO(app)
% Read data from Table
Table = app.UITable.Data;
% assignin('base','Table',Table)
% ListInd = logical(Table(:,1));
StructurePara = Table(:,2:end);

O.Phi_D1   = StructurePara(1,1);
O.Psi_D1   = StructurePara(1,2);
O.Theta_D1 = StructurePara(1,3);

O.Phi_D2   = StructurePara(2,1);
O.Psi_D2   = StructurePara(2,2);
O.Theta_D2 = StructurePara(2,3);

O.Displacement = StructurePara(2,4:6);

% Read ther GUI inputs
O.Rot_X = app.EditField_RotX.Value;
O.Rot_Y = app.EditField_RotY.Value;
O.Rot_Z = app.EditField_RotZ.Value;

% For Local Mode frequencies ----------------------------------------------
O.NLFreq  = app.EditField_LocFreq.Value;
O.LFreq   = app.EditField_LabelFreq.Value;
O.Anharm  = app.EditField_Anharm.Value;
O.L_Index = str2num(app.EditField_LabelInd.Value);

% For molecule ploting ----------------------------------------------------
O.Plot_Atoms = app.CheckBox_Atom.Value;
O.Plot_Bonds = app.CheckBox_Bond.Value;
O.Plot_Axis  = app.CheckBox_Axis.Value;