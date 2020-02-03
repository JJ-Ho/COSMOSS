%% Load Strucuture and find coordinate of O
S = Data_Comb2.Structure;
O_Ind = strcmp(S.AtomName,'O');
O_XYZ = S.XYZ(O_Ind,:);

%% Set grid spacing to be 0.1
XY = O_XYZ(:,1:2);
ShiftD = floor(min(XY(:))).*-1;

O_XYZ_Grid = XY+ShiftD;
O_XYZ_Grid = round(O_XYZ_Grid);

%% create a sparse matrix
M_L = max(O_XYZ_Grid(:));
M = sparse(O_XYZ_Grid(:,1),O_XYZ_Grid(:,2),O_XYZ(:,3),M_L,M_L);

%%
[X,Y]   = meshgrid(1:M_L,1:M_L);
[Xq,Yq] = meshgrid(1:0.1:M_L,1:0.1:M_L);
Z = full(M);

Zero_Height = 0;
Z(eq(Z,0))  = Zero_Height;
Max_Height  = max(Z(:));

% conv2
% Lineshape = 'Lorentzian';
Lineshape = 'Gaussian';
LineWidth = 3;
[ConvL2,~] = Conv_LineShape(2,Lineshape,1:M_L,LineWidth);
CVL_Z      = real(conv2(Z,ConvL2.lnshpf_R,'same'));
CVL_Z = CVL_Z.*(Max_Height/max(CVL_Z(:)));

ZShift = 3;
Zq = interp2(X,Y,CVL_Z,Xq,Yq,'Cubic');

X_Plot = Yq-ShiftD;
Y_Plot = Xq-ShiftD;
Z_Plot = Zq+ZShift;

hF = figure;
hAx = axes('Parent',hF);
hold(hAx,'on')
    %mesh(hAx,X_Plot,Y_Plot,Z_Plot)
    surf(hAx,X_Plot,Y_Plot,Z_Plot,'LineStyle','none','FaceAlpha',1)
    S.Draw(hAx)
hold(hAx,'off')

colormap(hAx,'cool')
caxis(hAx,[Zero_Height+0.5,Max_Height]+ZShift)

camlight 
% camlight('headlight')
camlight(-209,40) 
lighting gouraud

view([-137,23])