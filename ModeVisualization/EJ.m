function [EJ_T,Phi,Theta] = EJ(P,A,N_Grid)

P1 = P(1);
P2 = P(2);
P3 = P(3);
A1 = A(1);
A2 = A(2);
A3 = A(3);

% generate R(Phi,Theta)
phi   = linspace(0,2*pi,N_Grid);
theta = linspace(-pi/2,pi/2,N_Grid);
[Phi,Theta] = meshgrid(phi,theta);
T = pi/2-Theta(:); % make the vector perpendicular to incident beam direction
P = Phi(:);

t2 = cos(0);
t3 = sin(P);
t4 = cos(P);
t5 = cos(T);
t6 = sin(0);
t7 = sin(T);
R = reshape([-t3.*t6+t2.*t4.*t5,t4.*t6+t2.*t3.*t5,-t2.*t7,-t2.*t3-t4.*t5.*t6,t2.*t4-t3.*t5.*t6,t6.*t7,t4.*t7,t3.*t7,t5],[],3,3);
RT = reshape(permute(R,[3,1,2]),3,[]);

%%
% Signal(3)
E3 = [cosd(P3), sind(P3)];
J3 = [-cosd(A3), 0, sind(A3); % negative sign is for reflective geometry
              0, 1,       0];
EJ3 = reshape(E3*J3*RT,[],3);          
% EJ3 = reshape((R*(E3*J3)')',3,[]);          

% Visible(2)
E2 = [cosd(P2), sind(P2)];
J2 = [cosd(A2), 0, sind(A2);
             0, 1,       0];
EJ2 = reshape(E2*J2*RT,[],3);
% EJ2 = reshape((R*(E2*J2)')',3,[]);

% IR(1)
E1 = [cosd(P1), sind(P1)];
J1 = [cosd(A1), 0, sind(A1);
             0, 1,       0];
EJ1 = reshape(E1*J1*RT,[],3);  
% EJ1 = reshape((R*(E1*J1)')',3,[]); 


[Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3);
EJ_T = EJ3(:,Ja(:)).*EJ2(:,Jb(:)).*EJ1(:,Jc(:)); % [N_Grid^2 * 27]