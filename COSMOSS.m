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

% check if path is added otherwise, initailize the path
if ~eq(exist('TwoDSFG_Main','file'),2)
    Initialization
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @COSMOSS_OpeningFcn, ...
                   'gui_OutputFcn',  @COSMOSS_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
               
%- Check if GUI Layout Tool box exist ------------------------------------
T = ver;
UseLayoutToolBox = any(strcmp(cellstr(char(T.Name)), 'GUI Layout Toolbox'));

if and(UseLayoutToolBox,~nargin)
    gui_State.gui_LayoutFcn = @GUI_COSMOSS_Base;
    CreatMainGUI = 1;
else
    CreatMainGUI = 0;
end
% ------------------------------------------------------------------------   

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

% and(UseLayoutToolBox,any(isempty(gui_State.gui_Callback)))

if or(nargout,CreatMainGUI)
    if nargin
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        hMain = gui_mainfcn(gui_State, varargin{:});   
        varargout{1} = hMain;
    end
else
    gui_mainfcn(gui_State, varargin{:});
end

%- Call createInterface to create GUI elements and update handles ---------
if and(UseLayoutToolBox,~nargin)
    GUI_Main = GUI_COSMOSS(hMain);
    handles.hMain    = hMain;
    handles.GUI_Main = GUI_Main; % export GUI handles to handles
    guidata(hMain, handles);
end
% ------------------------------------------------------------------------
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

StructModel    = get(handles.GUI_Main.StructListBox,'Value');
[hStructure,~,~] = StructureModel(StructModel,handles.hMain);

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

GUI_Inputs = ParseGUI_Main(handles);
OneDSFG    = OneDSFG_Main(handles.Structure,GUI_Inputs);
PlotOneDSFG(OneDSFG,GUI_Inputs);

% Update share data
handles.OneDSFG = OneDSFG;
guidata(hObject,handles);

function TwoDIR_Callback(hObject, eventdata, handles)
% hObject    handle to Plot2DIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Read GUI
GUI_Inputs = ParseGUI_Main(handles);

%% Calculate TwoD response
Structure = handles.Structure;

if eq(GUI_Inputs.Sampling,1)
    % Pre-allocate
    GridSize     = length(GUI_Inputs.FreqRange);

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
    
    StandardDiv = GUI_Inputs.FWHM/(2*sqrt(2*log(2)));
    P_FlucCorr  = GUI_Inputs.P_FlucCorr/100; % turn percentage to number within 0~1
    
    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);
    
    
    for i = 1:GUI_Inputs.Sample_Num
        
        TSTART(i) = tic;
        
        % Add diagonal disorder
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv.*randn(1,1)*ones(Num_Modes,1);
        else 
            Fluctuation = StandardDiv.*randn(Num_Modes,1); 
        end
        Structure.freq = Freq_Orig + Fluctuation;
       
        [Tmp_SG,Tmp_Res] = TwoDIR_Main(Structure,GUI_Inputs);

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
    [SpectraGrid,Response] = TwoDIR_Main(Structure,GUI_Inputs);

    H     = Response.H;
    Mu_Ex = Response.Mu.Trans_Ex;
end

%% Conv2D linshape and make figure

CVL = Conv2D(SpectraGrid,GUI_Inputs);

% Make figure
% Plot2DIR(CVL,FreqRange,H,Mu_Ex,Daig_Cut)
Plot2DIR(CVL,GUI_Inputs.FreqRange,H,Mu_Ex)

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
GUI_Inputs = ParseGUI_Main(handles);

%% Calculate TwoD response
Structure = handles.Structure;

if eq(GUI_Inputs.Sampling,1)
    % Pre-allocate
    GridSize     = length(GUI_Inputs.FreqRange);

    Rephasing    = zeros(GridSize);
    NonRephasing = zeros(GridSize);
    SpecAccuR1   = zeros(GridSize);
    SpecAccuR2   = zeros(GridSize);
    SpecAccuR3   = zeros(GridSize);
    SpecAccuNR1  = zeros(GridSize);
    SpecAccuNR2  = zeros(GridSize);
    SpecAccuNR3  = zeros(GridSize);
    
    StandardDiv = GUI_Inputs.FWHM/(2*sqrt(2*log(2)));
    P_FlucCorr  = GUI_Inputs.P_FlucCorr/100; % turn percentage to number within 0~1
    
    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);
    
    for i = 1:GUI_Inputs.Sample_Num
        
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
        
        [Tmp_SG,Tmp_Res] = TwoDSFG_Main(Structure,GUI_Inputs);
        
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
    [SpectraGrid,Response] = TwoDSFG_Main(Structure,GUI_Inputs);
end
%% Covolution
CVL = Conv2D(SpectraGrid,GUI_Inputs);
  
%% Plot 2DSFG spectra
f = figure;
set(f,'Unit','normalized') % use normalized scale
CVLRS = -1.*real(CVL.selected);

contour(GUI_Inputs.FreqRange,GUI_Inputs.FreqRange,CVLRS,GUI_Inputs.Num_Contour,'LineWidth',2)

% Normalization
% CVLRSN = CVLRS ./max(abs(CVLRS(:)));
% contour(GUI_Inputs.FreqRange,GUI_Inputs.FreqRange,CVLRSN,GUI_Inputs.Num_Contour,'LineWidth',2)

% shading flat
ax = get(f,'CurrentAxes');
set(ax,'DataAspectRatio',[1 1 1])
% Plot diagonal line
hold on; plot(GUI_Inputs.FreqRange,GUI_Inputs.FreqRange); hold off

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

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'H', handles)
disp('Updated handles exported!')