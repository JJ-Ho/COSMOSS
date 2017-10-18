% SFG spectra simualtion with instant update of molecular orientation 
clear m
close all

GUI_Inputs.debug = '1';
GUI_Inputs.FreqRange = 1600:1700;

OneDSFG = OneDSFG_Main(Structure_Data,GUI_Inputs);
OneDSFG.SpecType     = 'SFG';
OneDSFG.FilesName    = 'SFG AR test';

E = OneDSFG.E;
J = OneDSFG.J;
Beta = OneDSFG.MolFrame;
Ex_F1 = OneDSFG.H.Sort_Ex_F1;
FreqRange = OneDSFG.freq_OneD;

%% prep for communication 
t=1;
t_max = 100;
SampleRate = 20; % Hz
XArray = 1:t_max;
O = zeros(length(XArray),3);

hF = figure;
SY = 25;
hAx_Spec = subplot(1,2,1);
hAx_Spec.YLim = [-SY,SY];

L = 20;
hAx_Mol = subplot(1,2,2);
hAx_Mol.XLim = [-L,L];
hAx_Mol.YLim = [-L,L];
hAx_Mol.ZLim = [-L,L];
view(hAx_Mol,[131,27])
rotate3d(hAx_Mol,'on')

m = mobiledev;

m.SampleRate = SampleRate;
m.OrientationSensorEnabled = 1;
m.Logging = 1;

pause(1)

%% Aqurire Orientation and generate response 
while(t<t_max)
    TSTART1 = tic;
    pause(1/SampleRate)
    
    [Read,T_Read] = orientlog(m);
    O(t,:) = Read(end,:)./180.*pi;
    
    
    % Spectrum genration
    R3_Mol2Lab = R3_ZYZ_0(O(t,1),O(t,2),O(t,3));
    
    Response = E * J * R3_Mol2Lab * Beta;
    AccuGrid = Bin1D(Ex_F1,Response,FreqRange);
    
    OneDSFG.Response1D = AccuGrid;
    
    Plot1D(hAx_Spec,OneDSFG,GUI_Inputs);
    hAx_Spec.YLim = [-SY,SY];
    
    TIME1 = toc(TSTART1);
    
    % Draw molecule
    TSTART2 = tic;
    
    XYZ_0 = Structure.XYZ;
    XYZ_0 = bsxfun(@minus,XYZ_0,sum(XYZ_0,1)./size(XYZ_0,1));
    
    R1_Mol2Lab = R1_zyx_0(O(t,1),O(t,2),O(t,3));
    XYZ_1 = (R1_Mol2Lab*XYZ_0')';
    
    Carbon_Pos   = XYZ_1(Structure.AtomSerNo(:,1),:);
    Oxygen_Pos   = XYZ_1(Structure.AtomSerNo(:,2),:);
    Nitrogen_Pos = XYZ_1(Structure.AtomSerNo(:,3),:);
    
    cla(hAx_Mol)
    axes(hAx_Mol)
    
    hold on
    Conn = Connectivity(XYZ_1);
    gplot3(Conn,XYZ_1);
    
    plot3(Carbon_Pos(:,1)  ,Carbon_Pos(:,2)  ,Carbon_Pos(:,3)  ,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerSize',5)
    plot3(Oxygen_Pos(:,1)  ,Oxygen_Pos(:,2)  ,Oxygen_Pos(:,3)  ,'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',5)
    plot3(Nitrogen_Pos(:,1),Nitrogen_Pos(:,2),Nitrogen_Pos(:,3),'LineStyle','none','Marker','o','MarkerFaceColor','b','MarkerSize',5)
    
    hold off

    
    hAx_Mol.XGrid ='on';
    hAx_Mol.YGrid ='on';
    hAx_Mol.ZGrid ='on';
    hAx_Mol.Box = 'on';
    axis equal
    hAx_Mol.XLim = [-L,L];
    hAx_Mol.YLim = [-L,L];
    hAx_Mol.ZLim = [-L,L];
    drawnow
    
    t = t +1;
    TIME2 = toc(TSTART2);
    disp([num2str(t) '-th loop finished...'])
    disp(['Calculation ' num2str(t) ' finished within '  num2str(TIME1) '...'])
    disp(['Drawing ' num2str(t) ' finished within '  num2str(TIME2) '!'])
end

m.Logging = 0;
m.OrientationSensorEnabled = 0;
clear m;

