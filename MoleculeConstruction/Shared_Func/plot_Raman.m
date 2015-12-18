function plot_Raman(Raman,Center,RT_scale,N_mesh)
% draw Raman tensor

% parameters
%RT_scale = 0.5;
%N_mesh = 20;

Orig_Frame = [1,0,0;0,1,0;0,0,1]';
    
% Extract ellipsoid info from Raman tensor
[V,D] = eig(Raman);
SemiAxisL = RT_scale.*diag(D);
E_Axis    = V;

% generate ellipsoid coordinate in Original Frame 
[Ex0, Ey0, Ez0] = ellipsoid(0,0,0,SemiAxisL(1),SemiAxisL(2),SemiAxisL(3),N_mesh);
OrigE = [Ex0(:),Ey0(:),Ez0(:)]; 

% rotation
RM = Euler_Rot(E_Axis,Orig_Frame);
Rot_XYZ = (RM*OrigE')';
Rot_XM = reshape(Rot_XYZ(:,1),N_mesh+1,N_mesh+1);
Rot_YM = reshape(Rot_XYZ(:,2),N_mesh+1,N_mesh+1);
Rot_ZM = reshape(Rot_XYZ(:,3),N_mesh+1,N_mesh+1);

% translation 
Trans_XM = Rot_XM + Center(1);
Trans_YM = Rot_YM + Center(2);
Trans_ZM = Rot_ZM + Center(3);

surf(Trans_XM,Trans_YM,Trans_ZM,...
    'LineStyle','-',...
    'EdgeAlpha',0.2,...
    'FaceColor',[204,152,255]./255,...
    'FaceAlpha',0.2);
