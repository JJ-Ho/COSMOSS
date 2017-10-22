function plot_Raman(hAx,Raman,Center,N_mesh,F_Color,Raman_Type)
%% debug
% clear all
% hF =figure;
% hAx = axes;
% Raman = [1,0,0;0,-5,0;0,0,10];
% Raman = [5.5,0,-4.5;0,-5,0;-4.5,0,5.5];
% Center = [0,0,0];
% RT_scale = 1;
% N_mesh = 20;
% F_Color = [1,0,0];
% Raman_Type = 1;

%% draw Raman tensor    
% Extract ellipsoid info from Raman tensor
[V,D] = eig(Raman);
% SemiAxisL = RT_scale.*diag(D);
SemiAxisL = diag(D);
RM = V; % Rotational matrix

switch Raman_Type
    case 'Arrow'
        %% Plot three 3D arrows with color representing sign
        Axes_Scale = abs(SemiAxisL)';
        Principle_Axes = bsxfun(@times,V,Axes_Scale);
        
        % Color of arrow
        ColorOrder = zeros(6,3);
        for i = 1:3
            if sign(SemiAxisL(i)) > 0
                ColorOrder(i,:) = [1,0,0];
                ColorOrder(i+3,:) = [1,0,0];
            else
                ColorOrder(i,:) = [0,0,1];
                ColorOrder(i+3,:) = [0,0,1];
            end
        end
        
        % Arror head
        Axes_Scale_N = Axes_Scale./max(Axes_Scale);
        Arror_W = 2*sqrt(Axes_Scale_N);
        Arror_H = 2*Arror_W;
        
        Arrow_End = [Principle_Axes,-Principle_Axes]'; % double side arrow
        Arrow_End = bsxfun(@plus,Arrow_End,Center); % move the ned point with Center translation
        Arror_Orig = bsxfun(@times,Center,ones(6,1));
        
        set(hAx,'ColorOrder',ColorOrder)
        arrow3(Arror_Orig,Arrow_End,'o2',Arror_W,Arror_H);
        
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
    
    SemiAxisL = abs(SemiAxisL);
    [D_x, D_y, D_z] = cylinder(Disk_Thickness_Array,N_mesh);  
    % XY disk
    Dxy = [   SemiAxisL(1).*D_x(:),   SemiAxisL(2).*D_y(:), Disk_Thickness.*D_z(:)];
    % YZ disk
    Dyz = [ Disk_Thickness.*D_z(:),   SemiAxisL(2).*D_y(:),   SemiAxisL(3).*D_x(:)];
    %Dyz = [ Disk_Thickness.*D_y(:),   SemiAxisL(2).*D_z(:),   SemiAxisL(3).*D_x(:)];
    % XY disk
    Dzx = [   SemiAxisL(1).*D_x(:), Disk_Thickness.*D_z(:),   SemiAxisL(3).*D_y(:)];
    
    hold on
    % rotation / translation and plot the disks
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
    patch(Trans_Dxy_XM',Trans_Dxy_YM',Trans_Dxy_ZM',F_Color(1,:),...
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
    patch(Trans_Dyz_XM',Trans_Dyz_YM',Trans_Dyz_ZM',F_Color(2,:),...
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
    patch(Trans_Dzx_XM',Trans_Dzx_YM',Trans_Dzx_ZM',F_Color(3,:),...
          'FaceLighting','gouraud',...
          'FaceAlpha',PatchFaceAlpha)
    hold off
    
    case 'Ellipsoid'
        %% Plot single ellipsoid
        % generate ellipsoid coordinate in Original Frame 
        [Ex0, Ey0, Ez0] = ellipsoid(0,0,0,SemiAxisL(1),SemiAxisL(2),SemiAxisL(3),N_mesh);
        OrigE = [Ex0(:),Ey0(:),Ez0(:)]; 

        % rotation
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
     
    case 'Hyper Ellipsoid'
        %% Plot hyper ellipsoid so that the radius = E'*Alpha*E
        phi   = linspace(0,2*pi,N_mesh);
        theta = linspace(-pi/2,pi/2,N_mesh);
        [Phi,Theta] = meshgrid(phi,theta);

        T = Theta(:);
        P = Phi(:);
        V1 = [cos(T).*cos(P),cos(T).*sin(P),sin(T)];
        V2 = V1;
        %V2 = [-sin(T).*cos(P),-sin(T).*sin(P),cos(T)]; % for cross polarization
        
        Rho = sum((V1*Raman).*V2,2);
        Rho = reshape(Rho,size(Theta));

        [X,Y,Z] = sph2cart(Phi,Theta,abs(Rho));
        X = X + Center(1);
        Y = Y + Center(2);
        Z = Z + Center(3);
        colormap('cool')
        
        caxis([-1,1])
        hSurf = surf(X,Y,Z,sign(Rho)); % colormapping sign only
        
        Transparency = 0.3;
        hSurf.EdgeColor = 'interp';
        hSurf.FaceAlpha = Transparency;
        hSurf.EdgeAlpha = Transparency;
end

%% debug
% axis equal
% view([40,20]) 
% lightangle(-45,30)
% hAx.XLabel.String = 'x';
% hAx.YLabel.String = 'y';
% hAx.ZLabel.String = 'z';
% hAx.Box = 'on';
% hAx.XLim = [-10,10];
% hAx.YLim = [-5,5];
% hAx.ZLim = [-10,10];
