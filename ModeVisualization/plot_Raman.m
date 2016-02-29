function plot_Raman(hAx,Raman,Center,RT_scale,N_mesh,F_Color)
%% debug
% hF =figure;
% hAx = axes;
% % Raman = [1,0,0;0,-5,0;0,0,10];
% Raman = [5.5,0,-4.5;0,-5,0;-4.5,0,5.5];
% Center = [0,0,0];
% RT_scale = 1;
% N_mesh = 20;
% F_Color = [1,0,0];

%% draw Raman tensor

% Draw_Type = 'Ball';
Draw_Type = 'Disk';

Orig_Frame = [1,0,0;0,1,0;0,0,1]';
    
% Extract ellipsoid info from Raman tensor
[V,D] = eig(Raman);
SemiAxisL = RT_scale.*diag(D);
E_Axis    = V;

% hold on
switch Draw_Type
    case 'Ball'
    %% Plot single ellipsoid
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

    surf(hAx,Trans_XM,Trans_YM,Trans_ZM,...
        'LineStyle','-',...
        'EdgeAlpha',0.2,...
        'FaceColor',F_Color,...
        'FaceAlpha',0.2);

    case 'Disk'
    %% Plot three ellipse disk
    Disk_Thickness = 0.05;
    Disk_Thickness_Array = ones(2,1);
    
    PatchFaceAlpha = 1.0;
    F_Color     = bsxfun(@times,ones(3),[1,0,1]);
    Cir_F_Color = bsxfun(@times,ones(3),[1,0,1]);
    for j = 1:3
        if eq(sign(SemiAxisL(j)),-1)
                F_Color(j,:) = [0,1,1];
            Cir_F_Color(j,:) = [0,1,1];
        end
    end
    
    [D_x, D_y, D_z] = cylinder(Disk_Thickness_Array,N_mesh);  
    % XY disk
    Dxy = [   SemiAxisL(1).*D_x(:),   SemiAxisL(2).*D_y(:), Disk_Thickness.*D_z(:)];
    % XY disk
    Dyz = [ Disk_Thickness.*D_z(:),   SemiAxisL(2).*D_y(:),   SemiAxisL(3).*D_x(:)];
    % XY disk
    Dzx = [   SemiAxisL(1).*D_x(:), Disk_Thickness.*D_z(:),   SemiAxisL(3).*D_y(:)];
    
    % rotation / translation and plot the disks
    RM = Euler_Rot(E_Axis,Orig_Frame);
    Rot_Dxy = (RM*Dxy')';
    Rot_Dxy_XM = reshape(Rot_Dxy(:,1),2,N_mesh+1);
    Rot_Dxy_YM = reshape(Rot_Dxy(:,2),2,N_mesh+1);
    Rot_Dxy_ZM = reshape(Rot_Dxy(:,3),2,N_mesh+1);
    Trans_Dxy_XM = Rot_Dxy_XM + Center(1);
    Trans_Dxy_YM = Rot_Dxy_YM + Center(2);
    Trans_Dxy_ZM = Rot_Dxy_ZM + Center(3);
    surf(hAx,Trans_Dxy_XM,Trans_Dxy_YM,Trans_Dxy_ZM,...
        'LineStyle','none',...
        'EdgeAlpha',0.2,...
        'FaceColor',Cir_F_Color(1,:),...
        'FaceAlpha',1);
    patch(Trans_Dxy_XM',Trans_Dxy_YM',Trans_Dxy_ZM',F_Color(3,:),...
          'FaceLighting','gouraud',...
          'FaceAlpha',PatchFaceAlpha)
    
    
    Rot_Dyz = (RM*Dyz')';
    Rot_Dyz_XM = reshape(Rot_Dyz(:,1),2,N_mesh+1);
    Rot_Dyz_YM = reshape(Rot_Dyz(:,2),2,N_mesh+1);
    Rot_Dyz_ZM = reshape(Rot_Dyz(:,3),2,N_mesh+1);
    Trans_Dyz_XM = Rot_Dyz_XM + Center(1);
    Trans_Dyz_YM = Rot_Dyz_YM + Center(2);
    Trans_Dyz_ZM = Rot_Dyz_ZM + Center(3);
    surf(hAx,Trans_Dyz_XM,Trans_Dyz_YM,Trans_Dyz_ZM,...
        'LineStyle','none',...
        'EdgeAlpha',0.2,...
        'FaceColor',Cir_F_Color(2,:),...
        'FaceAlpha',1);
    patch(Trans_Dyz_XM',Trans_Dyz_YM',Trans_Dyz_ZM',F_Color(1,:),...
          'FaceLighting','gouraud',...
          'FaceAlpha',PatchFaceAlpha)
    
    Rot_Dzx = (RM*Dzx')';
    Rot_Dzx_XM = reshape(Rot_Dzx(:,1),2,N_mesh+1);
    Rot_Dzx_YM = reshape(Rot_Dzx(:,2),2,N_mesh+1);
    Rot_Dzx_ZM = reshape(Rot_Dzx(:,3),2,N_mesh+1);
    Trans_Dzx_XM = Rot_Dzx_XM + Center(1);
    Trans_Dzx_YM = Rot_Dzx_YM + Center(2);
    Trans_Dzx_ZM = Rot_Dzx_ZM + Center(3);
    surf(hAx,Trans_Dzx_XM,Trans_Dzx_YM,Trans_Dzx_ZM,...
        'LineStyle','none',...
        'EdgeAlpha',0.2,...
        'FaceColor',Cir_F_Color(3,:),...
        'FaceAlpha',1);
    patch(Trans_Dzx_XM',Trans_Dzx_YM',Trans_Dzx_ZM',F_Color(2,:),...
          'FaceLighting','gouraud',...
          'FaceAlpha',PatchFaceAlpha)
    
    
end
% hold off

%% debug
% axis equal
% lightangle(-45,30)
