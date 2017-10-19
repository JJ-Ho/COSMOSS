function O = Parse_PDB_AmideOne(app)

O.Preprocessed = app.CheckBox_MDsnapshots.Value;
% Orientation -------------------------------------------------------------
O.Phi_D   =  app.EditField_Phi.Value;
O.Psi_D   =  app.EditField_Psi.Value;
O.Theta_D =  app.EditField_Theta.Value;

% For Local Mode frequencies ----------------------------------------------
O.NLFreq  = app.EditField_LocFreq.Value;
O.LFreq   = app.EditField_LabelFreq.Value;
O.Anharm  = app.EditField_Anharm.Value;
O.L_Index = str2num(app.EditField_LabelInd.Value);

% For molecule ploting ----------------------------------------------------
O.Plot_Atoms      = app.CheckBox_Atom.Value;
O.Plot_Bonds      = app.CheckBox_Bond.Value;
O.Plot_Axis       = app.CheckBox_Axis.Value;
O.Plot_Label      = app.CheckBox_Label.Value;
O.Plot_SideChain  = app.CheckBox_SideChain.Value;
