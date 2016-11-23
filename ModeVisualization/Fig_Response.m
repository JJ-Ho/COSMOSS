function Fig_Response(hAx, GUI_Inputs, Structure, OneDSFG, GUI_Data_hMain) 
% Plot hyper ellipsoid so that the radius = (ExJ)x(LxRxbeta)

%% define constants
N_Grid = 30;
ScaleFactor = 10;

A1 = GUI_Data_hMain.A_IR;
A2 = GUI_Data_hMain.A_Vis1D;
A3 = GUI_Data_hMain.A_Sig1D;
Ang = [A1,A2,A3];

P1 = GUI_Data_hMain.P_IR;
P2 = GUI_Data_hMain.P_Vis1D;
P3 = GUI_Data_hMain.P_Sig1D;
Polar = [P1,P2,P3];


hold on
%% Generate ExJ(psi,theta,Ai...,Pi...)
% The relative orientation of laser beams are determined by the incidnet 
% angles(Ai) and the polarization angle (Pi) defined on the incident beam 
% path. The molecular response to this laser setup configuration is a 
% function of the molecular orientation relative to the fixed experiment
% frame. Instead of viewing the response as function of molecule rotatein
% the exp. frame, we can fix the molecule and rotate the whole laser setup
% instead. These two picture are mathematically equivenlent but cause less
% computation if we apply the rotation on ExJ on each laer beam follow by 
% tensor product each of ExJ. This is effectly a [1x2]*[2x3]*[3x3]
% operation and follow by tensor product of however many [1x3] vectors. It
% is way faster than generating [3^Nx3^N] rotational matrix to rotate the
% higher order (N =3 for SFG, =5 for 2DSFG) molecular response. 
% If we further assume that the molecule system is azimutal symetric then 
% we can eliminate the phi angle out of three Euler angle and thus we can
% map the response (R) as function of two Euler angle so we can visulize it
% in (R,psi,theta) 3D coordinate. 

% % generate R(Psi,Theta)
% phi   = linspace(0,2*pi,N_Grid);
% theta = linspace(-pi/2,pi/2,N_Grid);
% [Phi,Theta] = meshgrid(phi,theta);
% 
% T = Theta(:);
% P = Phi(:);
% 
% 
% 
% V1 = [cos(T).*cos(P),cos(T).*sin(P),sin(T)];
% % V1 = [sin(T).*cos(P),sin(T).*sin(P),cos(T)];
% V2 = V1;
% V3 = V1;
% %V2 = [-sin(T).*cos(P),-sin(T).*sin(P),cos(T)]; % for cross polarization
% 
% [Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3);
% V3 = V3(:,Ja(:)).*V2(:,Jb(:)).*V1(:,Jc(:));

[EJ_M,Phi,Theta] = EJ(Polar,Ang,N_Grid);

%% Deal with response L x <R> x beta
% selecte mode
EigneVec_Ind = GUI_Inputs.EigneVec_Ind;

% Center of Exciton mode
EigVecM      = OneDSFG.H.Sort_Ex_V(2:end,2:end); % get ride of ground state
EigVecM2     = EigVecM.^2;
Center_Ex_MF = EigVecM2*(Structure.center);
Center       = Center_Ex_MF(EigneVec_Ind,:);

% Response 
% Response = OneDSFG.MolFrame;
Response = OneDSFG.LabFrame;
Rho = EJ_M*Response(:,EigneVec_Ind);
Rho = reshape(Rho,N_Grid,N_Grid);

% scale
Rho = Rho./ScaleFactor;

[X,Y,Z] = sph2cart(Phi,Theta,abs(Rho));
X = X + Center(1);
Y = Y + Center(2);
Z = Z + Center(3);
colormap('cool')

caxis([-1,1])
hSurf = surf(hAx,X,Y,Z,sign(Rho)); % colormapping sign only

Transparency = 0.5;
hSurf.EdgeColor = 'interp';
hSurf.FaceAlpha = Transparency;
hSurf.EdgeAlpha = Transparency;
hold off