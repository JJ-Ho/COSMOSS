function S_G09_LF = R_MF2LF(S_G09_MF,GUI_Inputs)
% This function update the new molecular frame and lab frame definition
% from GUI_Inputs

%% Process Inputs
Center_Ind = GUI_Inputs.MF_Center;
Z_i_Ind    = GUI_Inputs.MF_Zi;
Z_f_Ind    = GUI_Inputs.MF_Zf;
XZ_i_Ind   = GUI_Inputs.MF_XZi;
XZ_f_Ind   = GUI_Inputs.MF_XZf;
Frame_Type = GUI_Inputs.Frame_Type;
LF_Phi     = GUI_Inputs.LF_Phi./180*pi;
LF_Psi     = GUI_Inputs.LF_Psi./180*pi;
LF_Theta   = GUI_Inputs.LF_Theta./180*pi;

%% If the default molecular frame definition is different from what it read
% from G09 output, Set the assigned center to zero and rotated the molecule
% with the New [Orientation] Frame (NF) info
Z_Ind  = [ Z_i_Ind, Z_f_Ind];
XZ_Ind = [XZ_i_Ind,XZ_f_Ind];
S_G09_NF = SD_SetFrame(S_G09_MF,Center_Ind,Z_Ind,XZ_Ind);

% Change molecule definition from X-Z plane to Y-Z plane, if any
switch Frame_Type
    case 1 %Frame_Type = 'XZ';
    case 2 %Frame_Type = 'YZ';
        R_XZ2YZ  = R1_ZYZ_0(-pi/2,0,0);
        S_G09_NF = SD_Rot(S_G09_NF,R_XZ2YZ);
end

%% Rotate the molecule from molecule frame to lab frame
Rot2LF   = R1_ZYZ_0(LF_Phi,LF_Psi,LF_Theta);
S_G09_LF = SD_Rot(S_G09_NF,Rot2LF);
