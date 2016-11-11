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
        %varargout{1} = hMain;
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

% UIWAIT makes COSMOSS wait for user response (see UIRESUME)
% uiwait(handles.Main);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = COSMOSS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function onListSelection(hObject, eventdata, handles)

StructModel          = get(handles.GUI_Main.StructListBox,'Value');
[fhModel,ModelList,~] = StructureModel(StructModel);

guidata_Model = feval(fhModel,handles.hMain);
disp(['COSMOSS using model ' ModelList{StructModel}])

%- Data push to sub GUI test ----------------
% guidata_Model.Test = 'Test';
% guidata(guidata_Model.hModel,guidata_Model)
%- Data push to sub GUI test ----------------

% handles.Structure.hModel = hModel;
guidata(hObject,handles)

function FTIR_Callback(hObject, eventdata, handles)
% hObject    handle to PlotFTIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_Inputs = ParseGUI_Main(handles);
Structure = handles.Structure;

hF = figure;
hAx = axes;

if eq(GUI_Inputs.Sampling,1)
    % Pre-allocate
    GridSize   = length(GUI_Inputs.FreqRange);
    Num_Modes  = Structure.Num_Modes;
    Freq_Orig  = Structure.freq;
    Response1D = zeros(GridSize,1);
    
    
    StandardDiv = GUI_Inputs.FWHM./(2*sqrt(2*log(2)));
    P_FlucCorr  = GUI_Inputs.P_FlucCorr/100; % turn percentage to number within 0~1
    
    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);
    
    for i = 1:GUI_Inputs.Sample_Num
        DynamicUpdate = handles.GUI_Main.DynamicUpdate.Value;
        UpdateStatus  = handles.GUI_Main.UpdateStatus.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end
        
        TSTART(i) = tic;
        % Add diagonal disorder
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv'.*(randn(1,1).*ones(Num_Modes,1));
        else 
            Fluctuation = StandardDiv'.*randn(Num_Modes,1); 
        end
        Structure.freq = Freq_Orig + Fluctuation;
        FTIR = FTIR_Main(Structure,GUI_Inputs);
        
        % recursive part
        Response1D = Response1D + FTIR.Response1D; % note freq is binned and sported, so direct addition work
        FTIR.Response1D = Response1D;
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
        
        while ~eq(DynamicUpdate,0)
            FTIR.FilesName = [Structure.FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title
            Plot1D(hAx,FTIR,GUI_Inputs);
            drawnow
            DynamicUpdate = 0;
        end
    end
    
        Total_TIME = sum(TIME);
        disp(['Total time: ' num2str(Total_TIME)])
        
else
    FTIR = FTIR_Main(Structure,GUI_Inputs);
    Plot1D(hAx,FTIR,GUI_Inputs);
end

%% Update FTIR data into guidata 
handles.FTIR       = FTIR;
guidata(hObject,handles)

function SFG_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSFG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_Inputs = ParseGUI_Main(handles);
Structure = handles.Structure;

hF = figure;
hAx = axes;

if eq(GUI_Inputs.Sampling,1)
    % Pre-allocate
    GridSize   = length(GUI_Inputs.FreqRange);
    Num_Modes  = Structure.Num_Modes;
    Freq_Orig  = Structure.freq;
    Response1D = zeros(GridSize,1);
    
    
    StandardDiv = GUI_Inputs.FWHM./(2*sqrt(2*log(2)));
    P_FlucCorr  = GUI_Inputs.P_FlucCorr/100; % turn percentage to number within 0~1
    
    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);
    
    for i = 1:GUI_Inputs.Sample_Num
        DynamicUpdate = handles.GUI_Main.DynamicUpdate.Value;
        UpdateStatus  = handles.GUI_Main.UpdateStatus.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end
        
        TSTART(i) = tic;
        % Add diagonal disorder
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv'.*(randn(1,1).*ones(Num_Modes,1));
        else 
            Fluctuation = StandardDiv'.*randn(Num_Modes,1); 
        end
        Structure.freq = Freq_Orig + Fluctuation;
        OneDSFG = OneDSFG_Main(Structure,GUI_Inputs);
        
        % recursive part
        Response1D = Response1D + OneDSFG.Response1D; % note freq is binned and sported, so direct addition work
        OneDSFG.Response1D = Response1D;
        
        TIME(i) = toc(TSTART(i));
        disp(['Run ' num2str(i) ' finished within '  num2str(TIME(i)) '...'])
        
        while ~eq(DynamicUpdate,0)
            OneDSFG.FilesName = [Structure.FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title
            Plot1D(hAx,OneDSFG,GUI_Inputs);
            drawnow
            DynamicUpdate = 0;
        end
    end
    
        Total_TIME = sum(TIME);
        disp(['Total time: ' num2str(Total_TIME)])
        
else
    OneDSFG = OneDSFG_Main(Structure,GUI_Inputs);
    Plot1D(hAx,OneDSFG,GUI_Inputs);
end

%% Update share data
handles.OneDSFG = OneDSFG;
guidata(hObject,handles);

function TwoDIR_Callback(hObject, eventdata, handles)
% hObject    handle to Plot2DIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Read GUI
GUI_Inputs = ParseGUI_Main(handles);

DynamicUpdate = handles.GUI_Main.DynamicUpdate.Value;
if DynamicUpdate
    hF = figure;
end

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
    
    StandardDiv = GUI_Inputs.FWHM./(2*sqrt(2*log(2)));
    P_FlucCorr  = GUI_Inputs.P_FlucCorr/100; % turn percentage to number within 0~1
    
    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);
    
    for i = 1:GUI_Inputs.Sample_Num

        DynamicUpdate = handles.GUI_Main.DynamicUpdate.Value;
        UpdateStatus  = handles.GUI_Main.UpdateStatus.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end
        
        TSTART(i) = tic;
        
        % Add diagonal disorder
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv'.*(randn(1,1).*ones(Num_Modes,1));
        else 
            Fluctuation = StandardDiv'.*randn(Num_Modes,1); 
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
        
        SpectraGrid.Rephasing    = Rephasing    ;
        SpectraGrid.NonRephasing = NonRephasing ;
        SpectraGrid.SpecAccuR1   = SpecAccuR1   ;
        SpectraGrid.SpecAccuR2   = SpecAccuR2   ;
        SpectraGrid.SpecAccuR3   = SpecAccuR3   ;
        SpectraGrid.SpecAccuNR1  = SpecAccuNR1  ;
        SpectraGrid.SpecAccuNR2  = SpecAccuNR2  ;
        SpectraGrid.SpecAccuNR3  = SpecAccuNR3  ;
        Response = Tmp_Res;
        
        while ~eq(DynamicUpdate,0)
            CVL = Conv2D(SpectraGrid,GUI_Inputs);
            CVL.FilesName = [Structure.FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title
            Plot2DIR(hF,CVL,GUI_Inputs);
            drawnow
            DynamicUpdate = 0;
        end
    end
   
    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])
    
else
    [SpectraGrid,Response] = TwoDIR_Main(Structure,GUI_Inputs);
end

%% Conv2D linshape and make figure
hF_final = figure;
CVL = Conv2D(SpectraGrid,GUI_Inputs);

CVL.FilesName = Structure.FilesName; % pass filesname for figure title
% Make figure
Plot2DIR(hF_final,CVL,GUI_Inputs);

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

DynamicUpdate = handles.GUI_Main.DynamicUpdate.Value;
if DynamicUpdate
    hF = figure;
end

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
    
    StandardDiv = GUI_Inputs.FWHM./(2*sqrt(2*log(2)));
    P_FlucCorr  = GUI_Inputs.P_FlucCorr/100; % turn percentage to number within 0~1
    
    TSTART = zeros(GUI_Inputs.Sample_Num,1,'uint64');
    TIME   = zeros(GUI_Inputs.Sample_Num,1);
    
    for i = 1:GUI_Inputs.Sample_Num
        
        DynamicUpdate = handles.GUI_Main.DynamicUpdate.Value;
        UpdateStatus  = handles.GUI_Main.UpdateStatus.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end
        
        TSTART(i) = tic;
        
        % Add diagonal disorder
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv'.*(randn(1,1).*ones(Num_Modes,1));
        else 
            Fluctuation = StandardDiv'.*randn(Num_Modes,1); 
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
        
        SpectraGrid.Rephasing    = Rephasing    ;
        SpectraGrid.NonRephasing = NonRephasing ;
        SpectraGrid.SpecAccuR1   = SpecAccuR1   ;
        SpectraGrid.SpecAccuR2   = SpecAccuR2   ;
        SpectraGrid.SpecAccuR3   = SpecAccuR3   ;
        SpectraGrid.SpecAccuNR1  = SpecAccuNR1  ;
        SpectraGrid.SpecAccuNR2  = SpecAccuNR2  ;
        SpectraGrid.SpecAccuNR3  = SpecAccuNR3  ;
        Response = Tmp_Res;
        
        while ~eq(DynamicUpdate,0)
            CVL = Conv2D(SpectraGrid,GUI_Inputs);
            CVL.FilesName = [Structure.FilesName,' ',num2str(i),'-th run...']; % pass filesname for figure title
            Plot2DSFG(hF,CVL,GUI_Inputs);
            drawnow
            DynamicUpdate = 0;
        end
        
    end

    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])
    
else
    [SpectraGrid,Response] = TwoDSFG_Main(Structure,GUI_Inputs);
end

%% Covolution and make figure
hF_final = figure;
CVL = Conv2D(SpectraGrid,GUI_Inputs);

CVL.FilesName = Structure.FilesName; % pass filesname for figure title
Plot2DSFG(hF_final,CVL,GUI_Inputs);

%% update TwoDSFG_Response into guidata
TwoDSFG             = Response;
TwoDSFG.SpectraGrid = SpectraGrid;
TwoDSFG.CVL         = CVL;

handles.TwoDSFG = TwoDSFG;
guidata(hObject,handles);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hMain', handles)
disp('Updated handles exported!')
