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
Conc_Scaling = GUI_Inputs.Conc_Scaling;
Trans_X      = GUI_Inputs.Trans_X;
Trans_Y      = GUI_Inputs.Trans_Y;
Trans_Z      = GUI_Inputs.Trans_Z;
Rot_Phi      = GUI_Inputs.Rot_Phi/180*pi;
Rot_Psi      = GUI_Inputs.Rot_Psi/180*pi;
Rot_Theta    = GUI_Inputs.Rot_Theta/180*pi;

TransV = [Trans_X,Trans_Y,Trans_Z];

%% Move the molecules before merge the two Structures
% Move both CoM to [0,0,0]
S1_0 = SD_Trans(S1,-S1.CoM);
S2_0 = SD_Trans(S2,-S2.CoM);

R = R1_ZYZ_0(Rot_Phi,Rot_Psi,Rot_Theta);
S2_0_R = SD_Rot(S2_0,R);
S2_T_R = SD_Trans(S2_0_R,TransV);

%% Pass all the movement to the second structure's children for drawing purpose
N_S = length(S2_0.Children);
for i = 1:N_S
    C       = S2.Children(i);
    C_0     = SD_Trans(C,-C.CoM);
    C_0_R   = SD_Rot(C_0,R);
    C_T_R   = SD_Trans(C_0_R,TransV);
    C_T_R_0 = SD_Trans(C_T_R,C.CoM);
    S2_T_R.Children(i) = C_T_R_0;
end

%% Scale the transitions of the second molecule
S2_T_R_S = SD_ScaleTransitions(S2_T_R,Conc_Scaling);

%% Merge the two and Output
SC = SD_Comb2(S1_0,S2_T_R_S);

SC.hPlotFunc  = @PlotComb2;
SC.GUI_Inputs = GUI_Inputs;

% % Export into Structure so it can be passsed around different GUIs
% Structure.StrucData1  = StrucData1_0;
% Structure.StrucData2  = StrucData2_RT;