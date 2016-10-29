% SFG spectra simualtion with instant update of molecular orientation 
clear m
close all

GUI_Inputs.debug = '1';
GUI_Inputs.FreqRange = 1600:1700;
% GUI_Inputs.CouplingType = 'NN_Mix_TDC';
GUI_Inputs.CouplingType = 'Cho_APB';
GUI_Inputs.Beta_NN = 0.8;

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

% Default values
defaultCouplingType = 'TDC';
defaultAvg_Phi      = 0;
defaultAvg_Theta    = 0;
defaultAvg_Psi      = 0;
defaultAvg_Rot      = 1;
defaultAvg_Mirror   = 1;
defaultA_IR         = 90;
defaultA_Vis        = 90;
defaultA_Sum        = 90;
defaultP_IR         = 90;
defaultP_Vis        = 90;
defaultP_Sum        = 90;
defaultBeta_NN      = 0.8;
defaultFreqRange    = 1650:1750;

% add Optional inputs / Parameters
addOptional(INPUT,'CouplingType',defaultCouplingType);
addOptional(INPUT,'Avg_Phi'     ,defaultAvg_Phi);
addOptional(INPUT,'Avg_Theta'   ,defaultAvg_Theta);
addOptional(INPUT,'Avg_Psi'     ,defaultAvg_Psi);
addOptional(INPUT,'Avg_Rot'     ,defaultAvg_Rot);
addOptional(INPUT,'Avg_Mirror'  ,defaultAvg_Mirror);
addOptional(INPUT,'A_IR'        ,defaultA_IR);
addOptional(INPUT,'A_Vis'       ,defaultA_Vis);
addOptional(INPUT,'A_Sum'       ,defaultP_Sum);
addOptional(INPUT,'P_IR'        ,defaultP_IR);
addOptional(INPUT,'P_Vis'       ,defaultP_Vis);
addOptional(INPUT,'P_Sum'       ,defaultA_Sum);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);
addOptional(INPUT,'FreqRange'   ,defaultFreqRange);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names

CouplingType = INPUT.Results.CouplingType;
Avg_Phi      = INPUT.Results.Avg_Phi;
Avg_Theta    = INPUT.Results.Avg_Theta;
Avg_Psi      = INPUT.Results.Avg_Psi;
Avg_Rot      = INPUT.Results.Avg_Rot;
Avg_Mirror   = INPUT.Results.Avg_Mirror;
A_IR         = INPUT.Results.A_IR;
A_Vis        = INPUT.Results.A_Vis;
A_Sum        = INPUT.Results.A_Sum;
P_IR         = INPUT.Results.P_IR;
P_Vis        = INPUT.Results.P_Vis;
P_Sum        = INPUT.Results.P_Sum;
Beta_NN      = INPUT.Results.Beta_NN;
FreqRange    = INPUT.Results.FreqRange;

%% Every thing that can be prep in molecular frame 

H = ExcitonH(Structure,'ExMode','OneEx','CouplingType',CouplingType,'Beta_NN',Beta_NN);

Ex_Freq = H.Sort_Ex_Freq;

% construct mu,alpha
Mu    = MuAlphaGen(Structure,H,'Mode','Mu');
Alpha = MuAlphaGen(Structure,H,'Mode','Alpha');

Mu_Ex    = Mu.Trans_Ex;
Alpha_Ex = Alpha.Trans_Ex;

% Generate Molecular frame SFG Responses
Num_Modes = Structure.Num_Modes;

ResMolFrame = zeros(Num_Modes,3^3);

for N = 1:Num_Modes
    ResMolFrame(N,:) = kron(squeeze(Alpha_Ex(N+1,1,:)),squeeze(Mu_Ex(1,N+1,:)));
end

% Orientation = Orientation/180*pi; % turn to radius unit
Avg_Phi_R   =   Avg_Phi/180*pi;
Avg_Psi_R   =   Avg_Psi/180*pi;
Avg_Theta_R = Avg_Theta/180*pi;

switch Avg_Rot
        
    case 1 %'Phi' C_Inf

        R_Avg = R3_ZYZ_1(Avg_Psi_R,Avg_Theta_R);
        
    case 2 %'Psi'

        R_Avg = R3_ZYZ_2(Avg_Phi_R,Avg_Theta_R);
        
    case 3 %'{Phi,Psi}'

        R_Avg = R3_ZYZ_12(Avg_Theta_R);
        
    case 4 %'Isotropic'
        
        R_Avg = R3_ZYZ_123;
    
    case 5 %'No Average'
       
        R_Avg = R3_ZYZ_0(Avg_Phi_R,Avg_Psi_R,Avg_Theta_R);
    
        
    otherwise
        disp('Avg_Angle is not support, dont know how to apply Rotational average...')
end

% Decide Mirror planes
switch Avg_Mirror
    
    case 1 % no mirror plane
        V = [1;1;1];
        
        Mirror_Mask = kron(kron(V,V),V);
        
    case 2 % sigma v, X=-X, Y=-Y
        V1 = [-1; 1;1];
        V2 = [ 1;-1;1];
        
        Sigma_X = kron(kron(V1,V1),V1);
        Sigma_Y = kron(kron(V2,V2),V2);
        Sigma_X(eq(Sigma_X,-1)) = 0;
        Sigma_Y(eq(Sigma_Y,-1)) = 0;
        
        Mirror_Mask = and(Sigma_X,Sigma_Y);
end

Beta_Mol = (bsxfun(@times,R_Avg*ResMolFrame',Mirror_Mask))';

% Jones Matrix convert XYZ to PS frame
% 
% JLabFrame = [freq, ppp, pps, psp, pss, spp, sps, ssp, sss]
% 
% Turn degrees into radius
A_IR  =  A_IR/180*pi;
A_Vis = A_Vis/180*pi;
A_Sum = A_Sum/180*pi;

J = JonesRef3(A_Sum,A_Vis,A_IR);

% E part, Plarization of each incident beams
E = EPolar3(P_Sum,P_Vis,P_IR);

%% Prep initial SFG_Data
SFG_Data.SpecType     = 'SFG';
SFG_Data.Response1D   = E*J*Beta_Mol';
SFG_Data.freq_OneD    = FreqRange;
SFG_Data.FilesName    = 'SFG AR test';
SFG_Data.CouplingType = CouplingType;

%% Aqurire Orientation and generate response 
hF = figure;
SY = 25;
hAx_Spec = subplot(1,2,1);
hAx_Spec.YLim = [-SY,SY];

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

%% Read the data from Arduino
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
    
    Response = E*J*R3_Mol2Lab*Beta_Mol';
    ResponseGrid = Bin1D(Ex_Freq,Response,FreqRange);
    
    SFG_Data.Response1D = ResponseGrid;
    
    Plot1D(hAx_Spec,SFG_Data,GUI_Inputs);
    hAx_Spec.YLim = [-SY,SY];
    
    
    
    % Draw molecule
    XYZ_0 = Structure.XYZ;
    XYZ_0 = bsxfun(@minus,XYZ_0,sum(XYZ_0,1)./size(XYZ_0,1));
    
    R1_Mol2Lab = R1_Quad(V_Quad(i,1),-V_Quad(i,2),-V_Quad(i,3),-V_Quad(i,4));
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
    
    disp([num2str(i) '-th loop finished...'])
end

%% close communication
fclose(s);
delete(s)
clear s;

