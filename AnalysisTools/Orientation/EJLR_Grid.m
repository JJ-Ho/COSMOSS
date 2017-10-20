function EJLR=EJLR_Grid(P,A,N,Geometry)
% This function calculate the E*J*R(phi,theta,psi) matrix for molecular
% response visulization in tensor form. 

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
psi   = linspace(0,2*pi,N*2);

[Phi,Theta,Psi] = meshgrid(phi,theta,psi); % size(Phi(:)) = N*2N*2N = NG
% [Phi,Theta,Psi] = ndgrid(phi,theta,psi); % size(Phi(:)) = 2N*N*2N = NG

R1  = R1_ZYZ_0_ND(Phi(:),Psi(:),Theta(:)); % [NG,R1,R2]
R1T = reshape(permute(R1,[2,1,3]),3,[]);   % [R1,NG*R2]

%% EJLR
EJLR = reshape(E1*J1*R1T,[],3); % [NG * R2]
