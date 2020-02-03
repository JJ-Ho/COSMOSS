function [M,Phi,Theta]=EJRR_2DSFG(GUI_Data_hMain,N_Grid)
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
% Avg_Rot = 5;
% 
% N_Grid = 10;

%% Read parameters from COSMOSS GUI
A_Pump1 = GUI_Data_hMain.A_Pump1/180*pi;
A_Pump2 = GUI_Data_hMain.A_Pump2/180*pi;
A_Probe = GUI_Data_hMain.A_Probe/180*pi;
A_Vis2D = GUI_Data_hMain.A_Vis2D/180*pi;
A_Sig2D = GUI_Data_hMain.A_Sig2D/180*pi;

P_Pump1 = GUI_Data_hMain.P_Pump1/180*pi;
P_Pump2 = GUI_Data_hMain.P_Pump2/180*pi;
P_Probe = GUI_Data_hMain.P_Probe/180*pi;
P_Vis2D = GUI_Data_hMain.P_Vis2D/180*pi;
P_Sig2D = GUI_Data_hMain.P_Sig2D/180*pi;

Avg_Rot = GUI_Data_hMain.Avg_Rot;

%% E
E5 = EPolar5(P_Sig2D,P_Vis2D,P_Probe,P_Pump2,P_Pump1);

%% J
J5 = JonesRef5(A_Sig2D,A_Vis2D,A_Probe,A_Pump2,A_Pump1);

%% <R>
Dimension = 5; % for SFG

switch Avg_Rot
        
    case 1 %'Phi' C_Inf
        R_Avg = LabFrameAvg('C4',Dimension);
        
    case 4 %'Isotropic'
        R_Avg = LabFrameAvg('Isotropic',Dimension);
        
    case 5 %'No Average'
        R_Avg = LabFrameAvg('C1',Dimension);
        
    otherwise
        disp('Avg_Angle is not support, dont know how to apply Rotational average...')
        return
end

% Decide Mirror planes
% switch Avg_Mirror
%     
%     case 1 % no mirror plane
%         V = [1;1;1];
%         
%         Mirror_Mask = kron(kron(V,V),V);
% 
%     case 2 % sigma v, X=-X, Y=-Y
%         V1 = [-1; 1;1];
%         V2 = [ 1;-1;1];
%         
%         Sigma_X = kron(kron(V1,V1),V1);
%         Sigma_Y = kron(kron(V2,V2),V2);
%         Sigma_X(eq(Sigma_X,-1)) = 0;
%         Sigma_Y(eq(Sigma_Y,-1)) = 0;
%         
%         Mirror_Mask = and(Sigma_X,Sigma_Y);
% end

%% R in 3D space
% generate R(Phi,Theta)
phi     = linspace(0,2*pi,N_Grid*2);
theta   = linspace(0,  pi,N_Grid  );
[Phi,Theta] = meshgrid(phi,theta  );

R5 = R5_ZYZ_0_ND(Phi(:),0,Theta(:)); % [G^2,243,243]
% RT = reshape(permute(R3,[2,3,1]),243,[]);
RT = reshape(permute(R5,[3,2,1]),243,[]);

%% EJ<R>R
EJR_R = E5*J5*R_Avg*RT;
M = reshape(EJR_R,243,[])'; % [N_Grid^2 * 243]
