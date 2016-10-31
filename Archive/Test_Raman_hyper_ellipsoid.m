N_Grid = 50;

phi   = linspace(0,2*pi,N_Grid);
theta = linspace(-pi/2,pi/2,N_Grid);
[Phi,Theta] = meshgrid(phi,theta);

T = Theta(:);
P = Phi(:);
V1 = [cos(T).*cos(P),cos(T).*sin(P),sin(T)];
V2 = V1;
% V2 = [-sin(T).*cos(P),-sin(T).*sin(P),cos(T)];

% Raman = [1,0,0;0,-5,0;0,0,10];
Raman = [5.5,0,-4.5;0,-5,0;-4.5,0,5.5];

Rho = sum((V1*Raman).*V2,2);
Rho = reshape(Rho,size(Theta));

[X,Y,Z] = sph2cart(Phi,Theta,abs(Rho));

hF =figure;
hAx = axes;
% surf(X,Y,Z,sign(Z))
% surf(X,Y,Z,X.^2+Y.^2+Z.^2)
surf(X,Y,Z,Rho)

axis equal
view([40,20]) 
hAx.XLabel.String = 'x';
hAx.YLabel.String = 'y';
hAx.ZLabel.String = 'z';
hAx.Box = 'on';