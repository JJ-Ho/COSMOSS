function [SC,S1_New,S2_New] = Comb2(S1,S2,GUI_Inputs)
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

%% Move the molecules before merge the two Structure
% Move both COM to [0,0,0]
StrucData1_0 = SD_Trans(S1,-S1.COM);
StrucData2_0 = SD_Trans(S2,-S2.COM);

StrucData2_R  = SD_Rot(StrucData2_0,Rot_Phi,Rot_Psi,Rot_Theta);
StrucData2_RT = SD_Trans(StrucData2_R,TransV);

%% Merge the two and Output
S1_New = StrucData1_0;
S2_New = StrucData2_RT;
SC = SD_Comb2(S1_New,S2_New);

% % Export into Structure so it can be passsed around different GUIs
% Structure.StrucData1  = StrucData1_0;
% Structure.StrucData2  = StrucData2_RT;