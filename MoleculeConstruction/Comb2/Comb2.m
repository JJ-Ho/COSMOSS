function SC = Comb2(S1,S2,GUI_Inputs)
%% Debug
% hStruc1        = Data_Comb2.hStruc1;
% hStruc2        = Data_Comb2.hStruc2;
% StrucGUI_Data1 = guidata(hStruc1);
% StrucGUI_Data2 = guidata(hStruc2);
% StrucData1     = StrucGUI_Data1.Structure;
% StrucData2     = StrucGUI_Data2.Structure;
% GUI_Inputs = ParseGUI_Comb2(Data_Comb2.hGUIs);

%% Prep parameters
%Conc_Scaling = GUI_Inputs.Conc_Scaling;
Trans_X      = GUI_Inputs.Trans_X;
Trans_Y      = GUI_Inputs.Trans_Y;
Trans_Z      = GUI_Inputs.Trans_Z;
Rot_Phi      = GUI_Inputs.Rot_Phi/180*pi;
Rot_Psi      = GUI_Inputs.Rot_Psi/180*pi;
Rot_Theta    = GUI_Inputs.Rot_Theta/180*pi;

TransV = [Trans_X,Trans_Y,Trans_Z];

%% Move the molecules before merge the two Structures
% Move both CoM to [0,0,0]
S1 = SD_Trans(S1,-S1.CoM);
S2 = SD_Trans(S2,-S2.CoM);

S2 = SD_Rot(S2,Rot_Phi,Rot_Psi,Rot_Theta);
S2 = SD_Trans(S2,TransV);

%% Merge the two and Output
SC = SD_Comb2(S1,S2);

% % Export into Structure so it can be passsed around different GUIs
% Structure.StrucData1  = StrucData1_0;
% Structure.StrucData2  = StrucData2_RT;