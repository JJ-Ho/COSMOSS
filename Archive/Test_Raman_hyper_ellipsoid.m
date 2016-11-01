clear all

Mu1 = [0,0,1];
Mu2 = [0,1,0];
Mu3 = [0,1,0];

% Raman = [1,0,0;0,-5,0;0,0,10];
Raman = [5.5,0,-4.5;0,-5,0;-4.5,0,5.5];

N_Grid = 50;

phi   = linspace(0,2*pi,N_Grid);
theta = linspace(-pi/2,pi/2,N_Grid);
[Phi,Theta] = meshgrid(phi,theta);

T = Theta(:);
P = Phi(:);
V1 = [cos(T).*cos(P),cos(T).*sin(P),sin(T)];
V2 = V1;
% V2 = [-sin(T).*cos(P),-sin(T).*sin(P),cos(T)];

%% Tensor version, this can extend to N-dimesion 

% Raman only
% Expand index for tensor product between V1 and V2
% [I,J] = ndgrid(1:3,1:3);
% V1_V2_T = V1(:,I(:)).*V2(:,J(:)); % [N_Grid^2,3^2]
% 
% Rho = V1_V2_T*Raman(:); % [[N_Grid^2,3^2]] * [3^2,1]
% Rho = reshape(Rho,size(Theta)); % [N_Grid,N_Grid]

% Rmana*Mu 
% V3 = V2;
% [I_V1,I_V2,I_V3] = ndgrid(1:3,1:3,1:3);
% V1_V2_T = V1(:,I_V1(:)).*V2(:,I_V2(:)).*V3(:,I_V3(:)); % [N_Grid^2,3^3]
% 
% [I_R,I_M] = ndgrid(1:9,1:3);
% Chi = Raman(I_R(:)) .* Mu1(I_M(:))';
% 
% Rho = V1_V2_T*Chi; % [[N_Grid^2,3^3]] * [3^3,1]
% Rho = reshape(Rho,size(Theta)); % [N_Grid,N_Grid]

% Rmana*Mu*Mu*Mu
V3 = V2;
V4 = V2;
V5 = V2;

[I_V1,I_V2,I_V3,I_V4,I_V5] = ndgrid(1:3,1:3,1:3,1:3,1:3);
V1_V5_T = V1(:,I_V1(:)).*V2(:,I_V2(:)).*V3(:,I_V3(:)).*V4(:,I_V4(:)).*V5(:,I_V5(:)); % [N_Grid^2,3^5]

[I_R1,I_M2,I_M3,I_M4] = ndgrid(1:9,1:3,1:3,1:3);
Chi = Raman(I_R1(:)) .* Mu1(I_M2(:))' .* Mu2(I_M3(:))' .* Mu3(I_M4(:))';

Rho = V1_V5_T*Chi; % [[N_Grid^2,3^2]] * [3^2,1]
Rho = reshape(Rho,size(Theta)); % [N_Grid,N_Grid]


%% Original version, V'*Alpha*V
% Rho = sum((V1*Raman).*V2,2);
% Rho = reshape(Rho,size(Theta));

%% Making figure
[X,Y,Z] = sph2cart(Phi,Theta,abs(Rho));

hF =figure;
hAx = axes;

colormap('cool')
caxis([-1,1])

hSurf = surf(X,Y,Z,sign(Rho));

Transparency = 0.8;
% hSurf.EdgeColor = 'interp';
hSurf.FaceAlpha = Transparency;
hSurf.EdgeAlpha = Transparency;

axis equal
view([40,20]) 
hAx.XLabel.String = 'x';
hAx.YLabel.String = 'y';
hAx.ZLabel.String = 'z';
hAx.Box = 'on';