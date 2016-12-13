% function [M,Phi,Theta]=EJRR(GUI_Data_hMain,N_Grid)
%% Debug
A_IR    = 90./180*pi;
A_Vis1D = 90./180*pi;
A_Sig1D = 90./180*pi;

P_IR    = 0;
P_Vis1D = 0;
P_Sig1D = 0;

Avg_Rot = 5;

N_Grid = 10;

%% Read parameters from COSMOSS GUI
% A_IR    = GUI_Data_hMain.A_IR/180*pi;
% A_Vis1D = GUI_Data_hMain.A_Vis1D/180*pi;
% A_Sig1D = GUI_Data_hMain.A_Sig1D/180*pi;
% 
% P_IR    = GUI_Data_hMain.P_IR/180*pi;
% P_Vis1D = GUI_Data_hMain.P_Vis1D/180*pi;
% P_Sig1D = GUI_Data_hMain.P_Sig1D/180*pi;
% 
% Avg_Rot = GUI_Data_hMain.Avg_Rot;

%% E
E3 = EPolar3(P_Sig1D,P_Vis1D,P_IR);

%% J
J3 = JonesRef3(A_Sig1D,A_Vis1D,A_IR);

%% <R>
Dimension = 3; % for SFG

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
phi   = linspace(0,2*pi,N_Grid);
theta = linspace(-pi/2,pi/2,N_Grid);
[Phi,Theta] = meshgrid(phi,theta);
% T = pi/2-Theta(:); % make the vector perpendicular to incident beam direction
T = Theta(:);
P = Phi(:);

R3 = R3_ZYZ_0_ND(P,0,T); % [G^2,27,27]
RT = reshape(permute(R3,[2,3,1]),27,[]);
% RT = reshape(permute(R3,[2,1,3]),[],27);

%% EJ<R>R
EJR_R = E3*J3*R_Avg*RT;
M = reshape(EJR_R,27,[])'; % [N_Grid^2 * 27]
