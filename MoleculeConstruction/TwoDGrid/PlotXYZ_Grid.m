function hF = PlotXYZ_Grid(Structure)

%% Decide connectivity
Num_Modes = Structure.Num_Modes;
Atom_Num  = Structure.Monomer.Atom_Num;
XYZ       = Structure.XYZ;
Center    = Structure.center;

hF = figure;
hold on
    %% draw bonds
    Conn      = Connectivity(XYZ);
    % Add S-C connection for Ester
    for k = 1:Num_Modes
        Conn(13+(k-1)*Atom_Num,3+(k-1)*Atom_Num) = 1; 
    end
    gplot3(Conn,XYZ);

    %% Draw atoms
    Carbon_Pos   = XYZ(Structure.AtomSerNo(:,1),:);
    O_Double_Pos = XYZ(Structure.AtomSerNo(:,2),:);
    O_Single_Pos = XYZ(Structure.AtomSerNo(:,3),:);

    plot3(Center(:,1)       ,Center(:,2)       ,Center(:,3)      ,'LineStyle','none','Marker','d','MarkerFaceColor','w')
    plot3(Carbon_Pos(:,1)   ,Carbon_Pos(:,2)   ,Carbon_Pos(:,3)  ,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerSize',10)
    plot3(O_Double_Pos(:,1) ,O_Double_Pos(:,2) ,O_Double_Pos(:,3),'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',10)
    plot3(O_Single_Pos(:,1) ,O_Single_Pos(:,2) ,O_Single_Pos(:,3),'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',10)

    
%     %% draw transition dipoles
%     TDV_Scale = 0.1;
%     Mu_S = TDV_Scale .* Mu; % Scale TDV vector in plot
%     
%     quiver3(Center(:,1),Center(:,2),Center(:,3),...
%             Mu_S(:,1),Mu_S(:,2),Mu_S(:,3),0,...
%             'LineWidth',2,...
%             'Color',[255,128,0]./256);
%         
%     %% draw Raman tensor using plot_Raman
%     RT_scale = 0.5;
%     N_mesh   = 20;
% 
%     for i = 1: Num_Modes
%         Raman  = squeeze(RamanM(i,:,:));
%         plot_Raman(Raman,Center(i,:),RT_scale,N_mesh)
%     end

    hold off

%% Figure options
axis equal;
rotate3d on

xlabel('X')
ylabel('Y')
