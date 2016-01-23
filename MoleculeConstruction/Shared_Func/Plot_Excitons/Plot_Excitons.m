function varargout = Plot_Excitons(varargin)
% PLOT_EXCITONS MATLAB code for Plot_Excitons.fig
%      PLOT_EXCITONS, by itself, creates a new PLOT_EXCITONS or raises the existing
%      singleton*.
%
%      H = PLOT_EXCITONS returns the handle to a new PLOT_EXCITONS or the handle to
%      the existing singleton*.
%
%      PLOT_EXCITONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_EXCITONS.M with the given input arguments.
%
%      PLOT_EXCITONS('Property','Value',...) creates a new PLOT_EXCITONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Plot_Excitons_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Plot_Excitons_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Plot_Excitons

% Last Modified by GUIDE v2.5 22-Jan-2016 17:20:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Plot_Excitons_OpeningFcn, ...
                   'gui_OutputFcn',  @Plot_Excitons_OutputFcn, ...
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


% --- Executes just before Plot_Excitons is made visible.
function Plot_Excitons_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Plot_Excitons (see VARARGIN)

% Choose default command line output for Plot_Excitons
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Exciton= GUI_Plot_Excitons(hObject);
handles.GUI_Exciton = GUI_Exciton; % export GUI handles to handles

% Get Main function's handles
if nargin > 3    
    if ishandle(varargin{1}) 
       hStructure = varargin{1};
    end
else
    disp('Running in stand alone mode, using TCO mode for debug purpose')  
    hStructure = Model_TCO(handles);
end

% Retreive structure data from model GUI
Model_Data = guidata(hStructure);
Structure = Model_Data.Structure;

% Update handles structure
handles.hStructure = hStructure;
handles.Structure  = Structure;
guidata(hObject, handles);

% update exciton info
Update_Exciton(hObject, eventdata, handles)

% UIWAIT makes Plot_Excitons wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Plot_Excitons_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Update_Exciton(hObject, eventdata, handles)
% Run OneDSFG to get the corresponding mu and alpha of exciton modes
GUI_Inputs.test = 'test';
Structure = handles.Structure;
OneDSFG = OneDSFG_Main(Structure,GUI_Inputs);

Ex_Freq     = OneDSFG.H.Sort_Ex_Freq(2:end);
Ex_Mu       = squeeze(OneDSFG.Mu.Trans_Ex(1,2:end,:));
Ex_Mu_Int   = sqrt(sum(Ex_Mu.^2,2));
Ex_Alpha    = squeeze(OneDSFG.Alpha.Trans_Ex(1,2:end,:));
Ex_Alpha_ZZ = Ex_Alpha(:,9);

Mode_List = [Ex_Freq,Ex_Mu_Int,Ex_Alpha_ZZ];

% Update handles structure
handles.Mode_List = Mode_List;
guidata(hObject, handles);

% update the list on GUI
GUI_Exciton = handles.GUI_Exciton;
set(GUI_Exciton.ModeList,'Data',Mode_List)


function Update_Figure(hObject, eventdata, handles)

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'H', handles)
disp('Updated handles exported!')