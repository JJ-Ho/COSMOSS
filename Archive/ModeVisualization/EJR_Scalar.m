function [M2,M1,Phi,Theta]=EJR_Scalar(GUI_Data_hMain,N_Grid)
% This function calculate the E*J*R(phi,theta) matrix for Chi(2)
% visulization in tensor form. It is easy to apply ensemble average in
% tensor form so compare to function EJ, I include ensemble average here. 

%% Debug
% A_IR    = 90./180*pi;
% A_Vis1D = 90./180*pi;
% A_Sig1D = 90./180*pi;
% 
% P_IR    = 0;
% P_Vis1D = 0;
% P_Sig1D = 0;
% 
% N_Grid = 10;

%% Read parameters from COSMOSS GUI
A_IR    = GUI_Data_hMain.A_IR/180*pi;
A_Vis1D = GUI_Data_hMain.A_Vis1D/180*pi;
A_Sig1D = GUI_Data_hMain.A_Sig1D/180*pi;

P_IR    = GUI_Data_hMain.P_IR/180*pi;
P_Vis1D = GUI_Data_hMain.P_Vis1D/180*pi;
P_Sig1D = GUI_Data_hMain.P_Sig1D/180*pi;

%% E
E1 = EPolar1(P_IR);
E2 = EPolar2(P_Sig1D,P_Vis1D);

%% J
J1 = JonesRef1(A_IR);
J2 = JonesRef2(A_Sig1D,A_Vis1D);

%% R in 3D space
% generate R(Phi,Theta)
phi     = linspace(0,2*pi,N_Grid*2);
theta   = linspace(0,  pi,N_Grid  );
[Phi,Theta] = meshgrid(phi,theta  );

R1 = R1_ZYZ_0_ND(Phi(:),0,Theta(:)); % [N_Grid^2,3,3]
R2 = R2_ZYZ_0_ND(Phi(:),0,Theta(:)); % [N_Grid^2,9,9]

R1T = reshape(permute(R1,[3,2,1]),3,[]);
R2T = reshape(permute(R2,[3,2,1]),9,[]);

%% EJR
EJR1 = E1*J1*R1T;
EJR2 = E2*J2*R2T;

M1 = reshape(EJR1,3,[])'; % [N_Grid^2 * 3]
M2 = reshape(EJR2,9,[])'; % [N_Grid^2 * 9]
