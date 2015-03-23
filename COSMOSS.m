function varargout = COSMOSS(varargin)
% COSMOSS MATLAB code for COSMOSS.fig
%      COSMOSS, by itself, creates a new COSMOSS or raises the existing
%      singleton*.
%
%      H = COSMOSS returns the handle to a new COSMOSS or the handle to
%      the existing singleton*.
%
%      COSMOSS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COSMOSS.M with the given input arguments.
%
%      COSMOSS('Property','Value',...) creates a new COSMOSS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before COSMOSS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to COSMOSS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help COSMOSS

% Last Modified by GUIDE v2.5 01-Oct-2014 15:09:19

% initailize the path
Initialization

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @COSMOSS_OpeningFcn, ...
                   'gui_OutputFcn',  @COSMOSS_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before COSMOSS is made visible.
function COSMOSS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to COSMOSS (see VARARGIN)

% Choose default command line output for COSMOSS
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Main = GUI_COSMOSS(hObject);
handles.GUI_Main = GUI_Main; % export GUI handles to handles

% Test of create figure by GUI_COSMOSS(hObject);
% GUI = guihandles(hObject);
% handles.GUI = GUI; % export GUI handles to handles

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes COSMOSS wait for user response (see UIRESUME)
% uiwait(handles.Main);


% --- Outputs from this function are returned to the command line.
function varargout = COSMOSS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function onListSelection(hObject, eventdata, handles)

StructModel = get(handles.GUI_Main.StructListBox,'Value');
switch StructModel
    case 1
        hStructure = Model_TCO(handles.hMain);
    case 2 
        hStructure = Model_PDB_AmideI(handles.hMain);
end

handles.Structure.hStructure = hStructure;
guidata(hObject,handles)

function FTIR_Callback(hObject, eventdata, handles)
% hObject    handle to PlotFTIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

GUI_Inputs = ParseGUI_Main(handles);
FTIR       = PlotFTIR(handles.Structure,GUI_Inputs);

% Update share data 
handles.GUI_Inputs = GUI_Inputs;
handles.FTIR       = FTIR;
guidata(hObject,handles)

function SFG_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSFG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Read GUI
O = ParseGUI_Main(handles);

%% Main
OneDSFG = OneDSFG_Main(handles.Structure,...
                      'Label_Index' ,O.Label_Index,...
                      'Label_Freq'  ,O.Label_Freq,...
                      'Coupling'    ,O.Coupling,...
                      'Beta_NN'     ,O.Beta_NN,...
                      'Avg_Option'  ,O.Avg_Rot,...
                      'Mirror_Plane',O.Avg_Mirror,...
                      'A_IR'        ,O.A_IR,...
                      'A_Vis'       ,O.A_Vis1D,...
                      'A_Sum'       ,O.A_Sig1D,...
                      'P_IR'        ,O.P_IR,...
                      'P_Vis'       ,O.P_Vis1D,...
                      'P_Sum'       ,O.P_Sig1D);

PlotOneDSFG(OneDSFG,...
            'PlotStick'  ,O.PlotStick,...
            'F_Min'      ,O.F_Min,...
            'F_Max'      ,O.F_Max,...
            'LineWidth'  ,O.LineWidth,...
            'Signal_Type',O.Signal_Type,...
            'LineShape'  ,O.LineShape)

% Update share data
handles.OneDSFG = OneDSFG;
guidata(hObject,handles);

function TwoDIR_Callback(hObject, eventdata, handles)
% hObject    handle to Plot2DIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Read GUI
O = ParseGUI_Main(handles);

%% Calculate TwoD response
Structure = handles.Structure;
PolarAng  = [O.P_Pump1,O.P_Pump2,O.P_Probe,O.P_Sig2D];

if eq(O.Sampling,1)
    % Pre-allocate
    GridSize     = length(O.FreqRange);

    Rephasing    = zeros(GridSize);
    NonRephasing = zeros(GridSize);
    SpecAccuR1   = zeros(GridSize);
    SpecAccuR2   = zeros(GridSize);
    SpecAccuR3   = zeros(GridSize);
    SpecAccuNR1  = zeros(GridSize);
    SpecAccuNR2  = zeros(GridSize);
    SpecAccuNR3  = zeros(GridSize);
    
    Num_Modes = Structure.Num_Modes;
    Freq_Orig = Structure.freq;
    
    StandardDiv = O.FWHM/2*sqrt(2*log(2));
    P_FlucCorr  = O.P_FlucCorr/100; % turn percentage to number within 0~1
    
    TSTART = zeros(O.Sample_Num,1,'uint64');
    TIME   = zeros(O.Sample_Num,1);
    
    
    for i = 1:O.Sample_Num
        
        TSTART(i) = tic;
        
        % Add diagonal disorder
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv.*randn(1,1)*ones(Num_Modes,1);
        else 
            Fluctuation = StandardDiv.*randn(Num_Modes,1); 
        end
        Structure.freq = Freq_Orig + Fluctuation;
       
        [Tmp_SG,Tmp_Res] = TwoDIR_Main(Structure,...
                                       'FreqRange'  ,O.FreqRange,...
                                       'Label_Index',O.Label_Index,...
                                       'Label_Freq' ,O.Label_Freq,...
                                       'Coupling'   ,O.Coupling,...
                                       'Beta_NN'    ,O.Beta_NN,...
                                       'PolarAng'   ,  PolarAng);
        
        Rephasing    = Rephasing    + Tmp_SG.Rephasing   ;
        NonRephasing = NonRephasing + Tmp_SG.NonRephasing;
        SpecAccuR1   = SpecAccuR1   + Tmp_SG.SpecAccuR1  ;
        SpecAccuR2   = SpecAccuR2   + Tmp_SG.SpecAccuR2  ;
        SpecAccuR3   = SpecAccuR3   + Tmp_SG.SpecAccuR3  ;
        SpecAccuNR1  = SpecAccuNR1  + Tmp_SG.SpecAccuNR1 ;
        SpecAccuNR2  = SpecAccuNR2  + Tmp_SG.SpecAccuNR2 ;
        SpecAccuNR3  = SpecAccuNR3  + Tmp_SG.SpecAccuNR3 ;   
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
    end
    
    SpectraGrid.Rephasing    = Rephasing    ;
    SpectraGrid.NonRephasing = NonRephasing ;
    SpectraGrid.SpecAccuR1   = SpecAccuR1   ;
    SpectraGrid.SpecAccuR2   = SpecAccuR2   ;
    SpectraGrid.SpecAccuR3   = SpecAccuR3   ;
    SpectraGrid.SpecAccuNR1  = SpecAccuNR1  ;
    SpectraGrid.SpecAccuNR2  = SpecAccuNR2  ;
    SpectraGrid.SpecAccuNR3  = SpecAccuNR3  ;
    
    Response = Tmp_Res;
    
    H     = Tmp_Res.H;
    Mu_Ex = Tmp_Res.Mu.Trans_Ex;
   
    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])
    
else
    [SpectraGrid,Response] = TwoDIR_Main(Structure,...
                                        'FreqRange'  ,O.FreqRange,...
                                        'Label_Index',O.Label_Index,...
                                        'Label_Freq' ,O.Label_Freq,...
                                        'Coupling'   ,O.Coupling,...
                                        'Beta_NN'    ,O.Beta_NN,...
                                        'PolarAng'   ,  PolarAng);
    H     = Response.H;
    Mu_Ex = Response.Mu.Trans_Ex;
end

%% Conv2D linshape and make figure
CVL = Conv2D(SpectraGrid,...
            'FreqRange',O.FreqRange,...
            'LineShape',O.LineShape,...
            'LineWidth',O.LineWidth,...
            'SpecType' ,O.SpecType,...
            'Pathway'  ,O.Pathway);   

% Make figure
% Plot2DIR(CVL,FreqRange,H,Mu_Ex,Daig_Cut)
Plot2DIR(CVL,O.FreqRange,H,Mu_Ex)

%% update TwoDIR_Response into guidata
TwoDIR = Response;
TwoDIR.SpectraGrid = SpectraGrid;
TwoDIR.CVL         = CVL;

handles.TwoDIR = TwoDIR;
guidata(hObject,handles);

function TwoDSFG_Callback(hObject, eventdata, handles)
% hObject    handle to Plot2DSFG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Read GUI
O = ParseGUI_Main(handles);

%% Calculate TwoD response
Structure = handles.Structure;
PolarAng  = [O.P_Pump1,O.P_Pump2,O.P_Probe,O.P_Vis2D,O.P_Sig2D];
INCAng    = [O.A_Pump           ,O.A_Probe,O.A_Vis2D,O.A_Sig2D];

if eq(O.Sampling,1)
    % Pre-allocate
    GridSize     = length(O.FreqRange);

    Rephasing    = zeros(GridSize);
    NonRephasing = zeros(GridSize);
    SpecAccuR1   = zeros(GridSize);
    SpecAccuR2   = zeros(GridSize);
    SpecAccuR3   = zeros(GridSize);
    SpecAccuNR1  = zeros(GridSize);
    SpecAccuNR2  = zeros(GridSize);
    SpecAccuNR3  = zeros(GridSize);
    
    StandardDiv = O.FWHM/2*sqrt(2*log(2));
    P_FlucCorr  = O.P_FlucCorr/100; % turn percentage to number within 0~1
    
    TSTART = zeros(O.Sample_Num,1,'uint64');
    TIME   = zeros(O.Sample_Num,1);
    
    for i = 1:O.Sample_Num
        
        TSTART(i) = tic;
        
        % Sample structure with assigned fluctuation
        Num_Modes = Structure.Num_Modes;
        Freq_Orig = Structure.freq;
        
        % Add diagonal disorder
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv.*randn(1,1)*ones(Num_Modes,1);
        else 
            Fluctuation = StandardDiv.*randn(Num_Modes,1); 
        end
        Structure.freq = Freq_Orig + Fluctuation;
        
        [Tmp_SG,Tmp_Res] = TwoDSFG_AmideI_Main(Structure,...
                                               'FreqRange'  ,O.FreqRange,...
                                               'Label_Index',O.Label_Index,...
                                               'Label_Freq' ,O.Label_Freq,...
                                               'Coupling'   ,O.Coupling,...
                                               'Beta_NN'    ,O.Beta_NN,...
                                               'PolarAng'   ,  PolarAng,...
                                               'INCAng'     ,  INCAng,...
                                               'Avg_Rot'    ,O.Avg_Rot,...
                                               'Avg_Mirror' ,O.Avg_Mirror);
        
        Rephasing    = Rephasing    + Tmp_SG.Rephasing   ;
        NonRephasing = NonRephasing + Tmp_SG.NonRephasing;
        SpecAccuR1   = SpecAccuR1   + Tmp_SG.SpecAccuR1  ;
        SpecAccuR2   = SpecAccuR2   + Tmp_SG.SpecAccuR2  ;
        SpecAccuR3   = SpecAccuR3   + Tmp_SG.SpecAccuR3  ;
        SpecAccuNR1  = SpecAccuNR1  + Tmp_SG.SpecAccuNR1 ;
        SpecAccuNR2  = SpecAccuNR2  + Tmp_SG.SpecAccuNR2 ;
        SpecAccuNR3  = SpecAccuNR3  + Tmp_SG.SpecAccuNR3 ;   
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
    end
    
    SpectraGrid.Rephasing    = Rephasing    ;
    SpectraGrid.NonRephasing = NonRephasing ;
    SpectraGrid.SpecAccuR1   = SpecAccuR1   ;
    SpectraGrid.SpecAccuR2   = SpecAccuR2   ;
    SpectraGrid.SpecAccuR3   = SpecAccuR3   ;
    SpectraGrid.SpecAccuNR1  = SpecAccuNR1  ;
    SpectraGrid.SpecAccuNR2  = SpecAccuNR2  ;
    SpectraGrid.SpecAccuNR3  = SpecAccuNR3  ;
    
    Response = Tmp_Res;
    
    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])
    
else
    
    % Sample structure with assigned fluctuation
    [SpectraGrid,Response] = TwoDSFG_AmideI_Main(Structure,...
                                               'FreqRange'  ,O.FreqRange,...
                                               'Label_Index',O.Label_Index,...
                                               'Label_Freq' ,O.Label_Freq,...
                                               'Coupling'   ,O.Coupling,...
                                               'Beta_NN'    ,O.Beta_NN,...
                                               'PolarAng'   ,  PolarAng,...
                                               'INCAng'     ,  INCAng,...
                                               'Avg_Rot'    ,O.Avg_Rot,...
                                               'Avg_Mirror' ,O.Avg_Mirror);
end
%% Covolution
CVL = Conv2D(SpectraGrid,...
            'FreqRange',O.FreqRange,...
            'LineShape',O.LineShape,...
            'LineWidth',O.LineWidth,...
            'SpecType' ,O.SpecType,...
            'Pathway'  ,O.Pathway);
        
        
%% Plot 2DSFG spectra
f = figure;
set(f,'Unit','normalized') % use normalized scale
CVLRS = real(CVL.sum);
CVLRSN = CVLRS ./max(CVLRS(:));
contour(O.FreqRange,O.FreqRange,CVLRSN,O.Num_Contour,'LineWidth',2)
% shading flat
ax = get(f,'CurrentAxes');
set(ax,'DataAspectRatio',[1 1 1])
% Plot diagonal line
hold on; plot(O.FreqRange,O.FreqRange); hold off

% % Set colorbar
% MAP = custom_cmap(Num_Contour);
% colormap(MAP)

colorbar
Amp = max(abs(caxis));
caxis([-Amp Amp])

% Call pointer
S.fh = f;
S.ax = ax;
Pointer_N(S)

%% update TwoDSFG_Response into guidata
TwoDSFG = Response;
TwoDSFG.SpectraGrid = SpectraGrid;

handles.TwoDSFG = TwoDSFG;
guidata(hObject,handles);


