function EJLR=EJLR_Grid(P,A,N,Geometry)
% This function calculate the E*J*R(phi,theta) matrix for Chi(2)
% visulization in tensor form. It is easy to apply ensemble average in
% tensor form so compare to function EJ, I include ensemble average here. 

%% Debug
% % clear all
% P = 0;
% A = 0.5*pi;
% N = 10;
% Geometry = 'Reflective';

%% Main
E1 = EPolar1(P); % E

switch Geometry
    case 'Transimisive'
        J1 = JonesTrans1(A); % J
    case 'Reflective'
        J1 = JonesRef1(A); % J
end

% Generate R(Phi,Psi,Theta)
phi   = linspace(0,2*pi,N*2);
theta = linspace(0,  pi,N  );
% psi   = linspace(0,2*pi,N*2);
psi = 0;
[Phi,Theta,Psi] = meshgrid(phi,theta,psi);

R1 = R1_ZYZ_0_ND(Phi(:),Psi(:),Theta(:)); % [4*N_Grid^3,3,3]
R1T = reshape(permute(R1,[3,2,1]),3,[]);

%% EJLR
EJLR = reshape(E1*J1*R1T,3,[])'; % [4*N_Grid^3 * 3]

