function varargout = Plot_Modes(varargin)
% PLOT_MODES MATLAB code for Plot_Modes.fig
%      PLOT_MODES, by itself, creates a new PLOT_MODES or raises the existing
%      singleton*.
%
%      H = PLOT_MODES returns the handle to a new PLOT_MODES or the handle to
%      the existing singleton*.
%
%      PLOT_MODES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_MODES.M with the given input arguments.
%
%      PLOT_MODES('Property','Value',...) creates a new PLOT_MODES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Plot_Modes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Plot_Modes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Plot_Modes

% Last Modified by GUIDE v2.5 24-Jan-2016 09:53:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Plot_Modes_OpeningFcn, ...
                   'gui_OutputFcn',  @Plot_Modes_OutputFcn, ...
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


% --- Executes just before Plot_Modes is made visible.
function Plot_Modes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Plot_Modes (see VARARGIN)

% Choose default command line output for Plot_Modes
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Modes = GUI_Plot_Modes(hObject);
handles.GUI_Modes = GUI_Modes; % export GUI handles to handles

% Get Structural modeling GUI's handles
if nargin > 3    
    if ishandle(varargin{1}) 
       hModel = varargin{1};
    end
else
    disp('Running in stand alone mode, using TCO modle for debug purpose')  
    hModel = Model_TCO(handles);
end

% Change Names on Plot_Exciton GUI to identify which Structural model is
% using
Model_Name = hModel.Name;
Plot_Modes_GUI_Name = GUI_Modes.hPlot_Modes.Name;
GUI_Modes.hPlot_Modes.Name = [Plot_Modes_GUI_Name, ': ', Model_Name];

% Export the handle of Structure Modle to guidata of Plot_Exciton
handles.hModel = hModel;
guidata(hObject, handles);

% update exciton info in Plot_Exciton
Update_Modes(hObject, eventdata, handles)

% UIWAIT makes Plot_Modes wait for user response (see UIRESUME)
% uiwait(handles.hPlot_Modes);


% --- Outputs from this function are returned to the command line.
function varargout = Plot_Modes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Update_Modes(hObject, eventdata, handles)
% Run OneDSFG to get the corresponding mu and alpha of exciton modes
% Retrieve Label index and Coupling model from COSMOSS GUI if any, if
% running Plot_Exciton stand alone for debug testing, give a field 'debug'
% to use the default values in OneDSFG_Main.m
GUI_Data_hModel = guidata(handles.hModel);
if isfield(GUI_Data_hModel,'hMain')
    GUI_Data_hMain = guidata(GUI_Data_hModel.hMain);
    MainGUI_Inputs = ParseGUI_Main(GUI_Data_hMain);
    %disp('Using the labeing index and coupling info from Main GUI')
else
    MainGUI_Inputs.debug = 'debug';
    GUI_Data_hModel.hMain = 'debug';
    %disp('Labeling index and coupling info come from defulat setting of OneDSFG_Main.m')
end

Structure = GUI_Data_hModel.Structure;
OneDSFG = OneDSFG_Main(Structure,MainGUI_Inputs);

Ex_Freq     = OneDSFG.H.Sort_Ex_Freq(2:end);
Ex_Mu       = squeeze(OneDSFG.Mu.Trans_Ex(1,2:end,:));
Ex_Mu_Int   = sqrt(sum(Ex_Mu.^2,2));
Ex_Alpha    = squeeze(OneDSFG.Alpha.Trans_Ex(1,2:end,:));
Ex_Alpha_ZZ = Ex_Alpha(:,9);

Mode_List = [Ex_Freq,Ex_Mu_Int,Ex_Alpha_ZZ];

% Update handles structure
handles.hMain          = GUI_Data_hModel.hMain;
handles.Structure      = Structure;
handles.MianGUI_Inputs = MainGUI_Inputs;
handles.OneDSFG        = OneDSFG;
guidata(hObject, handles);

% update the list on hPlot_Exciton GUI
set(handles.GUI_Modes.ModeList,'Data',Mode_List)

function Update_Figure(hObject, eventdata, handles)
%% Re assign variable names of Inputs
GUI_handle = handles.hGUI_Plot_Modes.Ex_Mode_Ind;
Mode_Ind = str2num(GUI_handle.String)+1; % shift by 1 to avoid ground state 

Structure = handles.Structure;
OneDSFG   = handles.OneDSFG;

% Mu     = OneDSFG.Mu.Trans_Ex(Mode_Ind,:);
% Alpha  = OneDSFG.Alpha.Trans_Ex(Mode_Ind,:);
% 
% Center = Structure.center(Mode_Ind,:);

%% Calculate the exciton center

%% update handles
handles.Mode_Ind = Mode_Ind;

guidata(hObject,handles)

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hPlot_Modes', handles)
disp('Updated handles exported!')