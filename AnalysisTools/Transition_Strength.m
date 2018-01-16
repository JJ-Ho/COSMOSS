function Transition_Strength(SpecData)
% Transition_Strength takes Hamiltonian (H) and Transition vectors (TV) to
% draw the energy level diagram with different opacity. The darker a state
% is the larger the transition strength is. 

%% Debug
% SpecData.SpecType = 'SFG';
% SpecData.H = Data_COSMOSS.Data_SFG.H;
% SpecData.Alpha = Data_COSMOSS.Data_SFG.Alpha;

%% Determine spectral Type
switch SpecData.SpecType
    case 'FTIR'
        ExType = '1ex';
        H  = SpecData.H;
        TV = SpecData.Mu;
    case 'SFG'
        ExType = '1ex';
        H  = SpecData.H;
        TV = SpecData.Alpha;
    case '2DIR'
        ExType = '2ex';
        H  = SpecData.Response.H;
        TV = SpecData.Response.Mu;
    case '2DSFG'
        ExType = '2ex';
        H  = SpecData.Response.H;
        TV = SpecData.Response.Alpha;
end
%% Figure parameters
% Box dimension to represent transition strength
Box_X = 0.6;
Box1_Y = 5;
Box2_Y = 5;
Color1 = [1,0,0]; % One Exciton modes
Color2 = [0,0,1]; % Two Exciton modes
Y_Blank = 5; % for Y axis range

%% One Exciton
N1 = H.Nmodes;
F1 = H.Sort_Ex_F1;
M1 = TV.M_Ex_01;
switch size(M1,2)
    case 3
        TV_Type = 'IR';
    case 9
        TV_Type = 'Raman';
end

M1L = TV.M_Ex_01_N;
M1L_N = M1L./max(M1L);

X1_V = (1:N1)';
X1_L_L = X1_V - Box_X/2;
X1_R_L = X1_V + Box_X/2;
X1_M_L = [X1_L_L,X1_R_L]';

X1_L_P = X1_V - Box_X/2.*0.8;
X1_R_P = X1_V + Box_X/2.*0.8;
X1_M_P = [X1_L_P,X1_L_P,X1_R_P,X1_R_P];

Y1_V = F1;
Y1_T = Y1_V  + Box1_Y/2.*M1L_N;
Y1_B = Y1_V  - Box1_Y/2.*M1L_N;

Y1_M_L = [Y1_V,Y1_V]';
Y1_M_P = [Y1_B,Y1_T,Y1_T,Y1_B];

% Patch Color
C1 = bsxfun(@times,ones(N1,1),Color1);

%% Two Exciton
switch ExType
    case '1ex'
        hF = figure;

        % One Exciton part
        hAx1 = axes('Parent',hF);
        for i = 1:N1
            patch(X1_M_P(i,:),Y1_M_P(i,:),C1(i,:),'FaceAlpha',M1L_N(i),'EdgeColor','none')
        end
        line(hAx1,X1_M_L,Y1_M_L,'Color','k')
        
        hAx1.XLim = [0,N1+1];
        hAx1.YLim = [min(Y1_V)-Y_Blank,max(Y1_V)+Y_Blank];
        hAx1.XTick = 1:N1;
        hAx1.YGrid = 'on';
        hAx1.XLabel.String = '1-Ex Mode #';
        hAx1.YLabel.String = 'F_1 (cm^{-1})';
        hAx1.FontSize = 16;
        hAx1.Title.String = [TV_Type ' transition strength'];
        
    case '2ex'
        N2 = H.StatesNum - H.Nmodes -1;
        F2 = H.Sort_Ex_F2;
        
        M2L = TV.M_Ex_12_N;
        M2L = M2L(:);
        M2L_N = M2L./max(M2L);

        X2_M = bsxfun(@times,(1:N1)',ones(1,N2));
        X2_V = X2_M(:);

        X2_L_L1 = X2_V - Box_X/2.*0.7;
        X2_R_L1 = X2_V - Box_X/2;
        X2_M_L1 = [X2_L_L1,X2_R_L1]';

        X2_L_L2 = X2_V + Box_X/2.*0.7;
        X2_R_L2 = X2_V + Box_X/2;
        X2_M_L2 = [X2_L_L2,X2_R_L2]';

        X2_L_P = X2_V - Box_X/2.*0.8;
        X2_R_P = X2_V + Box_X/2.*0.8;
        X2_M_P = [X2_L_P,X2_L_P,X2_R_P,X2_R_P];

        Y2_M = bsxfun(@minus,F2',2.*F1);
        Y2_V = Y2_M(:);
        Y2_T = Y2_V - Box2_Y/2*M2L_N;
        Y2_B = Y2_V + Box2_Y/2*M2L_N;

        Y2_M_L = [Y2_V,Y2_V]';
        Y2_M_P = [Y2_B,Y2_T,Y2_T,Y2_B];

        % Patch Color
        C2 = bsxfun(@times,ones(N1*N2,1),Color2);

        % Prepare the connection lines between the same F2 states
        X2_Link = [reshape(X2_M(1:end-1,:)+ Box_X/2,1,[]);...
                   reshape(X2_M(2:end  ,:)- Box_X/2,1,[])];

        Y2_Link = [reshape(Y2_M(1:end-1,:),1,[]);...
                   reshape(Y2_M(2:end  ,:),1,[])];

        %% Draw figure
        hF = figure;

        % One Exciton part
        hAx1 = subplot(4,1,4,'Parent',hF);
        for i = 1:N1
            patch(X1_M_P(i,:),Y1_M_P(i,:),C1(i,:),'FaceAlpha',M1L_N(i),'EdgeColor','none')
        end
        line(hAx1,X1_M_L,Y1_M_L,'Color','k')

        % Two Exciton part
        hAx2 = subplot(4,1,1:3,'Parent',hF);
        for j = 1:N1*N2
            patch(X2_M_P(j,:),Y2_M_P(j,:),C2(j,:),'FaceAlpha',M2L_N(j),'EdgeColor','none')
        end
        line(hAx2,X2_M_L1,Y2_M_L,'Color','k')
        line(hAx2,X2_M_L2,Y2_M_L,'Color','k')

        % Draw F2 state link lines
        line(hAx2,X2_Link,Y2_Link,'Color','k','LineStyle',':')

        %% Figure adjustment 
        hAx1.XLim = [0,N1+1];
        hAx1.YLim = [min(Y1_V)-Y_Blank,max(Y1_V)+Y_Blank];
        hAx1.XTick = 1:N1;
        hAx1.YGrid = 'on';
        hAx1.XLabel.String = '1-Ex Mode #';
        hAx1.YLabel.String = 'F_1 (cm^{-1})';
        hAx1.FontSize = 16;

        % Figure title
        hAx2.Title.String = [TV_Type ' transition strength'];
        hAx2.XLim = [0,N1+1];
        hAx2.YLim = [min(Y2_V)-Y_Blank,max(Y2_V)+Y_Blank];
        hAx2.YGrid = 'on';
        hAx2.XTick = [];
        hAx2.YLabel.String = '2-Ex shift, F_2 - 2F_1 (cm^{-1})';
        hAx2.FontSize = 16;
end