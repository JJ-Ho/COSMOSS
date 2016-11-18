function Fig_Response(hAx, GUI_Inputs, Structure, OneDSFG)       
hold on
%% Plot hyper ellipsoid so that the radius = E'*Alpha*E
        N_Grid = 30;

        phi   = linspace(0,2*pi,N_Grid);
        theta = linspace(-pi/2,pi/2,N_Grid);
        [Phi,Theta] = meshgrid(phi,theta);

        T = Theta(:);
        P = Phi(:);
        V1 = [cos(T).*cos(P),cos(T).*sin(P),sin(T)];
        V2 = V1;
        V3 = V1;
        %V2 = [-sin(T).*cos(P),-sin(T).*sin(P),cos(T)]; % for cross polarization

        [Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3);
        V3 = V3(:,Ja(:)).*V2(:,Jb(:)).*V1(:,Jc(:));
     
        % selecte mode
        EigneVec_Ind = GUI_Inputs.EigneVec_Ind;
        Center       = Structure.center(EigneVec_Ind,:);
        
        Response = OneDSFG.MolFrame;
        Rho = V3*Response(:,EigneVec_Ind);
        Rho = reshape(Rho,size(Theta));
        
        % scale
        Rho = Rho;
        
        [X,Y,Z] = sph2cart(Phi,Theta,abs(Rho));
        X = X + Center(1);
        Y = Y + Center(2);
        Z = Z + Center(3);
        colormap('cool')
        
        caxis([-1,1])
        hSurf = surf(hAx,X,Y,Z,sign(Rho)); % colormapping sign only
        
        Transparency = 0.3;
        hSurf.EdgeColor = 'interp';
        hSurf.FaceAlpha = Transparency;
        hSurf.EdgeAlpha = Transparency;
hold off