% function MixSpec(hF1,hF2,Ratio)
%% Debug
hF1 = figure(1);
hF2 = figure(2);
hF3 = figure(10);
hF4 = figure(4);
%% Main
hC1 = findobj(hF1,'Type','Contour');
hC2 = findobj(hF2,'Type','Contour');
hC3 = findobj(hF3,'Type','Contour');
hC4 = findobj(hF4,'Type','Contour');

hC1 = hC1(1);
hC2 = hC2(1);
hC3 = hC3(1);
hC4 = hC4(1);

Z1 = hC1.ZData;
Z2 = hC2.ZData;
Z3 = hC3.ZData;
Z4 = hC4.ZData;

R1 = 1;
R2 = 0.3;
R3 = 1;
R4 = 0.2;

Z = R1.*Z1 + R2.*Z2 + R3.*Z3 + R4.*Z4;
% Z = R1.*Z1 + R2.*Z2 + R3.*Z3;
X = hC1.XData;
Y = hC1.YData;

Num_Contour = 40;

%% Make figure
hF  = figure;
hAx = axes('Parent',hF);
hold(hAx,'on')
        [~,hC] = contourf(hAx,X,Y,Z,Num_Contour,'LineWidth',2);
        
        hC.LineStyle = 'none';
        C_MAP = SelectColormap('JetWhiteWide');      
        C_Amp = max(abs(Z(:)));
        DiagColor = 'k';
        
        Fig_ContourLine = 1;
        if Fig_ContourLine
            % copy the contour and make it into line only
            hCL = copyobj(hC,hAx);
            hCL.Fill = 'off';
            hCL.LineStyle = '-';
            hCL.LineWidth = 0.5;
            hCL.LineColor = [0,0,0]; 
            NC_half = floor(Num_Contour/2);
            CLevelList = linspace(-1,1,NC_half)*C_Amp;
            %Fig_CL_D_Ind = [floor(NC_half/2):floor(NC_half/2)+1];
            Fig_CL_D_Ind = floor(NC_half/2)+1;
            CLevelList(Fig_CL_D_Ind)=[]; % Deal with deleting zero contour line
            hCL.LevelList = CLevelList;
        end
colorbar(hAx)
colormap(hAx,C_MAP)
caxis(hAx,[-C_Amp,C_Amp])

% Plot diagonal line
plot(hAx,X,Y,'Color',DiagColor,'LineStyle','--')
hold(hAx,'off')

hAx.DataAspectRatio = [1,1,1];
hAx.FontSize = 14;
hAx.XLabel.String = 'Probe (cm^{-1})';
hAx.YLabel.String = 'Pump (cm^{-1})';
hAx.XLim = [1550,1800];
hAx.YLim = [1550,1800];