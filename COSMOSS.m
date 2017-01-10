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
                   'gui_LayoutFcn',  @GUI_Base_COSMOSS, ...
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

function hCOSMOSS = GUI_Base_COSMOSS(Singleton)
% this function create base figure. To utilize the original GUI building
% mechanism of Matlab and avoid Matlab default way to export GUI element 
% handles, this function only create the base layer with no GUI elements.
% Aftter calling "gui_mainfcn" we will call "GUI_COSMOSS" to build GUI
% elements.

% Create base figure
hCOSMOSS = figure;

hCOSMOSS.Units            = 'Pixels';
hCOSMOSS.Position         = [2 406 560 600];
hCOSMOSS.Name             = 'COSMOSS';
hCOSMOSS.ToolBar          = 'none';
hCOSMOSS.MenuBar          = 'none';
hCOSMOSS.NumberTitle      = 'off';
hCOSMOSS.IntegerHandle    = 'off';
hCOSMOSS.Tag              = 'hMain'; % tag to distinguish type of GUI
hCOSMOSS.HandleVisibility = 'Callback';

gui_Options.syscolorfig = 1;
setappdata(hCOSMOSS,'GUIDEOptions',gui_Options);

disp('Creating COSMOSS GUI Using GUI Layout Toolbox!')
disp('...')

function COSMOSS_OpeningFcn(hCOSMOSS, eventdata, GUI_data, varargin)
% This function has no output args, see OutputFcn.
% hCOSMOSS   handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% GUIdata    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to COSMOSS (see VARARGIN)

%- Check if GUI Layout Tool box exist ------------------------------------
T = ver;
UseLayoutToolBox = any(strcmp(cellstr(char(T.Name)), 'GUI Layout Toolbox'));

if UseLayoutToolBox
    hGUIs = GUI_COSMOSS(hCOSMOSS);
else
    diap('Please install GUI Layout toolbox...')
    return
end
% ------------------------------------------------------------------------   

% Prep necessary data to be export
GUI_data.hCOSMOSS = hCOSMOSS; 
GUI_data.hGUIs    = hGUIs; % export GUI handles to handles

guidata(hCOSMOSS, GUI_data);

function varargout = COSMOSS_OutputFcn(hCOSMOSS, eventdata, GUI_data) 
% varargout  cell array for returning output args (see VARARGOUT);
% hCOSMOSS    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% GUIdata    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hCOSMOSS;

function Export_Handle_Callback(hObject, eventdata, GUI_data)
% export handles back to work space
assignin('base', 'Data_COSMOSS', GUI_data)
disp('Updated GUI Data_COSMOSS exported!')
%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


function onListSelection(hObject, eventdata, GUI_data)
StructModel           = get(GUI_data.hGUIs.StructListBox,'Value');
[fhModel,ModelList,~] = StructureModel(StructModel);

hModel = feval(fhModel,'COSMOSS',GUI_data.hCOSMOSS);
disp(['COSMOSS using model ' ModelList{StructModel}])

% Pass hCOSMOSS to the fig file of each sub-GUIs so they can always access
% COSMOSS if needed
hModel.UserData = GUI_data.hCOSMOSS;

GUI_data.hModel = hModel;
guidata(hObject,GUI_data)

function FTIR_Callback(hObject, eventdata, GUI_data)
% update the laser seeting tab
hGUIs = GUI_data.hGUIs;
hGUIs.LaserSetting.SelectedTab = hGUIs.Tab_1D;

%% Main
GUI_Inputs = ParseGUI_Main(GUI_data.hGUIs);
Structure  = GUI_data.Structure;

hF  = figure;
hAx = axes('Parent',hF);

if eq(GUI_Inputs.Sampling,1) % GUI: Diag. Disorder
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
        DynamicUpdate = GUI_data.hGUIs.DynamicUpdate.Value;
        UpdateStatus  = GUI_data.hGUIs.UpdateStatus.Value;
        if and(~eq(i,1), and(eq(DynamicUpdate,1),~eq(UpdateStatus,1)))
            break
        end
        
        TSTART(i) = tic;
        
        S = Structure;
        % Replace mu and center if there are mutiple MD sanpshots
        if isfield(S,'N_File') && gt(S.N_File,1)
            S.mu = squeeze(S.mu(:,:,i));
            S.center = S.center(:,:,i);
        end
        % Add diagonal disorder if any
        Correlation_Dice = rand;

        if Correlation_Dice < P_FlucCorr
            Fluctuation = StandardDiv'.*(randn(1,1).*ones(Num_Modes,1));
        else 
            Fluctuation = StandardDiv'.*randn(Num_Modes,1); 
        end
        S.freq = Freq_Orig + Fluctuation;
        
        % Calculate FTIR
        FTIR = FTIR_Main(S,GUI_Inputs);
        
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
    
    Plot1D(hAx,FTIR,GUI_Inputs);
    Total_TIME = sum(TIME);
    disp(['Total time: ' num2str(Total_TIME)])
        
else
    FTIR = FTIR_Main(Structure,GUI_Inputs);
    Plot1D(hAx,FTIR,GUI_Inputs);
end

%% Update FTIR data into guidata 
GUI_data.FTIR = FTIR;
guidata(hObject,GUI_data)

function SFG_Callback(hObject, eventdata, GUI_data)
% update the laser seeting tab
hGUIs = GUI_data.hGUIs;
hGUIs.LaserSetting.SelectedTab = hGUIs.Tab_1D;

%% Main
GUI_Inputs = ParseGUI_Main(GUI_data.hGUIs);
Structure  = GUI_data.Structure;

hF  = figure;
hAx = axes('Parent',hF);

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
        DynamicUpdate = GUI_data.hGUIs.DynamicUpdate.Value;
        UpdateStatus  = GUI_data.hGUIs.UpdateStatus.Value;
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
    
        Plot1D(hAx,OneDSFG,GUI_Inputs);
        Total_TIME = sum(TIME);
        disp(['Total time: ' num2str(Total_TIME)])
        
else
    OneDSFG = OneDSFG_Main(Structure,GUI_Inputs);
    Plot1D(hAx,OneDSFG,GUI_Inputs);
end

%% Update share data
GUI_data.OneDSFG = OneDSFG;
guidata(hObject,GUI_data);

function TwoDIR_Callback(hObject, eventdata, GUI_data)
%% Read GUI
hGUIs = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_Main(hGUIs);

% update the laser seeting tab
hGUIs.LaserSetting.SelectedTab = hGUIs.Tab_2D;

DynamicUpdate = GUI_Inputs.DynamicUpdate;
if DynamicUpdate
    hF = figure;
end

%% Calculate TwoD response
Structure = GUI_data.Structure;

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

        DynamicUpdate = hGUIs.DynamicUpdate.Value; % directly access the GUI elment so can get the most recnt values
        UpdateStatus  = hGUIs.UpdateStatus.Value;
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
       
        %[Tmp_SG,Tmp_Res] = TwoDIR_Main(Structure,GUI_Inputs);
        [Tmp_SG,Tmp_Res] = TwoDIR_Main_Sparse(Structure,GUI_Inputs);
        
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
    %[SpectraGrid,Response] = TwoDIR_Main(Structure,GUI_Inputs);
    [SpectraGrid,Response] = TwoDIR_Main_Sparse(Structure,GUI_Inputs);
end

%% Conv2D linshape and make figure
hF_final = figure;
CVL = Conv2D(SpectraGrid,GUI_Inputs);

CVL.FilesName = Structure.FilesName; % pass filesname for figure title
% Make figure
Plot2DIR(hF_final,CVL,GUI_Inputs);

%% update TwoDIR_Response into guidata
TwoDIR             = Response;
TwoDIR.SpectraGrid = SpectraGrid;
TwoDIR.CVL         = CVL;
TwoDIR.SpecType   = '2DIR';

GUI_data.TwoDIR = TwoDIR;
guidata(hObject,GUI_data);

function TwoDSFG_Callback(hObject, eventdata, GUI_data)
%% Read GUI
hGUIs = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_Main(hGUIs);

% update the laser seeting tab
hGUIs.LaserSetting.SelectedTab = hGUIs.Tab_2D;

DynamicUpdate = GUI_Inputs.DynamicUpdate;
if DynamicUpdate
    hF = figure;
end

%% Calculate TwoD response
Structure = GUI_data.Structure;

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
        
        DynamicUpdate = hGUIs.DynamicUpdate.Value;
        UpdateStatus  = hGUIs.UpdateStatus.Value;
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
        % disp(num2str(Freq_Orig + Fluctuation)) %[debug]
        
        %[Tmp_SG,Tmp_Res] = TwoDSFG_Main(Structure,GUI_Inputs);
        [Tmp_SG,Tmp_Res] = TwoDSFG_Main_Sparse(Structure,GUI_Inputs);
        
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
    %[SpectraGrid,Response] = TwoDSFG_Main(Structure,GUI_Inputs);
    [SpectraGrid,Response] = TwoDSFG_Main_Sparse(Structure,GUI_Inputs);
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
TwoDSFG.SpecType    = '2DSFG';

GUI_data.TwoDSFG = TwoDSFG;
guidata(hObject,GUI_data);



function HCut_Callback(hObject, eventdata, GUI_data)
GUI_Inputs = ParseGUI_Main(GUI_data.hGUIs);

SpecType = GUI_data.hGUIs.AnalysisTools.SelectedTab.Title;
switch SpecType
    case '2DIR'
        Cut_F     = GUI_Inputs.HCut_2DIR;
        TwoD_Data = GUI_data.TwoDIR;
    case '2DSFG'
        Cut_F     = GUI_Inputs.HCut_2DSFG;
        TwoD_Data = GUI_data.TwoDSFG;
    otherwise
        disp('Spectra Type not supported')
end

hF_HCut = H_Cut(Cut_F,GUI_Inputs,TwoD_Data);

HCut.hF_HCut = hF_HCut;

%% Update share data
GUI_data.HCut = HCut;
guidata(hObject,GUI_data);