function hF = PlotCavity(hAx,obj_SD)
%% Make mirror surface
R = 40;
r = 1:0.5:10;
phi = 0:1:360;
Rim = max(r);

[Radius,PHI] = meshgrid(r,phi);
X = R-sqrt(R^2 - Radius.^2);
Y = Radius.*cos(PHI);
Z = Radius.*sin(PHI);

dist = 20;
Mirror1_X =  X - dist;
Mirror2_X = -X + dist;

% optical axis
X_axis = -dist:1:dist;
Y_axis = zeros(size(X_axis));
Z_axis = zeros(size(X_axis));

% rim of the mirror
X_r = (R-sqrt(R^2 - Rim.^2)).*ones(size(phi));
Y_r = Rim.*cos(phi);
Z_r = Rim.*sin(phi);

%% Prep hAx
% hAx = 'new';
if ~ishandle(hAx)
    hF = figure; 
    hAx = axes('Parent',hF);
else
    hF = hAx.Parent;
end

hold(hAx,'on')
hM1 = surf(hAx,Mirror1_X,Y,Z);
hM2 = surf(hAx,Mirror2_X,Y,Z);

hOAxis = line(hAx,X_axis,Y_axis,Z_axis);
hRim1  = line(hAx, X_r - dist,Y_r,Z_r);
hRim2  = line(hAx,-X_r + dist,Y_r,Z_r);

hold(hAx,'off')

%% figure setting 
% Optical axis
hOAxis.LineWidth = 2;

% Rim
% hRim1.LineWidth = 2;
% hRim2.LineWidth = 2;

% mirrors
hM1.EdgeColor = 'interp';
hM2.EdgeColor = 'interp';
hM1.CData = X;
hM2.CData = X;
colormap(hAx,'copper')

% Axis
hAx.XLabel.String = 'X';
hAx.YLabel.String = 'Y';
hAx.ZLabel.String = 'Z';
hAx.XLim = [-dist,dist].*1.5;
hAx.YLim = [-Rim,Rim].*1.5;
hAx.ZLim = [-Rim,Rim].*1.5;
axis(hAx,'image');
rotate3d(hAx,'on')
view(hAx,[50,10])
axis(hAx,'image');
camlight
daspect([1 1 1]);
