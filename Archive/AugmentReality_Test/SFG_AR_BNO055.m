function SFG_AR_BNO055(Data_COSMOSS)
% SFG spectra simualtion with instant update of molecular orientation 

%% close communication and clean variables
if exist('s','var')
    fclose(s);
    delete(s)
end

%% 1DSFG 
G = ParseGUI_Main(Data_COSMOSS.hGUIs);
S = Data_COSMOSS.Structure;
SFG_Data = OneDSFG_Main(S,G);

E = SFG_Data.E; 
J = SFG_Data.Jones;
Beta_Mol = SFG_Data.LabFrame;

%% Aqurire Orientation and generate response 
hF = figure;
hAx_Spec = subplot(1,2,1);

L = 20;
hAx_Mol = subplot(1,2,2);
hAx_Mol.XLim = [-L,L];
hAx_Mol.YLim = [-L,L];
hAx_Mol.ZLim = [-L,L];
hAx_Mol.XLabel.String = 'X';
hAx_Mol.YLabel.String = 'Y';
hAx_Mol.ZLabel.String = 'Z';
view(hAx_Mol,[131,27])
rotate3d(hAx_Mol,'on')

%% Create serial object for Arduino
[~,Port] = unix('ls /dev/tty.usb*');
Port(Port==10) = []; % remove neline in unix, in Windows us "13"
% s = serial('/dev/tty.usbmodem1421'); % change the COM Port number as needed
% s = serial('/dev/tty.usbmodemFA131'); % change the COM Port number as needed
s = serial(Port);
s.BaudRate = 115200;
% Connect the serial port to Arduino
try
    fopen(s);
catch err
    fclose(instrfind);
    error('Make sure you select the correct COM Port where the Arduino is connected.');
end

%% handshaking between Arduino and Matlab
a = 'b';
while (a~='a') 
    a=fread(s,1,'uchar');
end
if (a=='a')
    disp('Serial read');
end

fprintf(s,'%c','a');
mbox = msgbox('Serial Communication setup'); uiwait(mbox);
fscanf(s,'%u');

%% Read the data from Arduino and calculate spectrum
t_max = 100;
V_Quad = zeros(1,4);

i = 1; % loop index
T_loop = 0;
T_loop_end = 0;

tic % Start timer
while(toc<t_max)

    i = i + 1;
    T_loop_end(i) = toc;
    T_loop(i) = toc - T_loop_end(i-1);
    
    try
        V_Quad(i,:) = readQuad(s,'R');
    catch
        disp('skip empty output...')
    end
    
    % normalize V_Quad
    V_Quad(i,:) = V_Quad(i,:)./norm(V_Quad(i,:));
    
    % Spectrum genration
    R3_Mol2Lab = R3_Quad(V_Quad(i,1),-V_Quad(i,2),-V_Quad(i,3),-V_Quad(i,4));
    
    Response = E*J*R3_Mol2Lab*Beta_Mol;
    ResponseGrid = Bin1D(SFG_Data.H.Sort_Ex_F1,Response,G.FreqRange);
    
    SFG_Data.Response1D = ResponseGrid;
    
    Plot1D(hAx_Spec,SFG_Data,G);
    SY = max(abs(ResponseGrid(:))) * 1.1 ;
    hAx_Spec.YLim = [-SY,SY];
    
    % Draw molecule
    XYZ_0 = S.XYZ;
    XYZ_0 = bsxfun(@minus,XYZ_0,sum(XYZ_0,1)./size(XYZ_0,1));
    
    R1_Mol2Lab = R1_Quad(V_Quad(i,1),-V_Quad(i,2),-V_Quad(i,3),-V_Quad(i,4));
    XYZ_1 = (R1_Mol2Lab*XYZ_0')';
    
    C_Ind = strcmp(S.AtomName,'C');
    O_Ind = strcmp(S.AtomName,'O');
    N_Ind = strcmp(S.AtomName,'N');
    
    Carbon_Pos   = XYZ_1(C_Ind,:);
    Oxygen_Pos   = XYZ_1(O_Ind,:);
    Nitrogen_Pos = XYZ_1(N_Ind,:);
    
    cla(hAx_Mol)
    axes(hAx_Mol)
    
    hold on
    Conn = Connectivity(S.AtomName,XYZ_1);
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
    
    disp([num2str(i) '-th loop finished...'])
end

%% close communication
fclose(s);
delete(s)
clear s;

