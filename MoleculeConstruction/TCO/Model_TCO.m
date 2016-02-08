function varargout = Model_TCO(varargin)
% Model_TCO MATLAB code for Model_TCO.fig
%      Model_TCO, by itself, creates a new Model_TCO or raises the existing
%      singleton*.
%
%      H = Model_TCO returns the handle to a new Model_TCO or the handle to
%      the existing singleton*.
%
%      Model_TCO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Model_TCO.M with the given input arguments.
%
%      Model_TCO('Property','Value',...) creates a new Model_TCO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_TCO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_TCO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Model_TCO

% Last Modified by GUIDE v2.5 01-Oct-2014 16:16:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Model_TCO_OpeningFcn, ...
                   'gui_OutputFcn',  @Model_TCO_OutputFcn, ...
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


% --- Executes just before Model_TCO is made visible.
function Model_TCO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_TCO (see VARARGIN)

% Choose default command line output for Model_TCO
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Struc = GUI_TCO(hObject);
handles.GUI_Struc = GUI_Struc; % export GUI handles to handles

% Get Main function's handles
if nargin > 3    
    if ishandle(varargin{1})
       hMain = varargin{1};
       Data_Main = guidata(hMain);

       handles.hMain = hMain;
       handles.Data_Main = Data_Main;
       
       % Reset Non-Label Frequency, anharmonicity, and F_min/F_Max to fit amideI mode
       GUI_Main  = Data_Main.GUI_Main;
       set(GUI_Main.NLFreq ,'String','1700')
       set(GUI_Main.LFreq  ,'String','1680')
       set(GUI_Main.Anharm ,'String','20')
       set(GUI_Main.Beta_NN,'String','5')
       set(GUI_Main.X_Min  ,'String','1650')
       set(GUI_Main.X_Max  ,'String','1750')
    end
else
    disp('Running Model_TCO in stand alone mode.')
end

% Update handles structure
guidata(hObject, handles);

% Generate the default structure
UpdateStructure(hObject, eventdata, handles)

% UIWAIT makes Model_TCO wait for user response (see UIRESUME)
% uiwait(handles.Model_TCO);


% --- Outputs from this function are returned to the command line.
function varargout = Model_TCO_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function UpdateStructure(hObject, eventdata, handles)
%% Read GUI variables
GUI_Struc = handles.GUI_Struc;
GUI_Inputs = ParseGUI_TCO(GUI_Struc);

%% Construct molecule
Structure = GetAcid(GUI_Inputs);

% Export into Structure so it can be passsed around different GUIs
Structure.StructModel = 1;                

%% Export result to Main guidata
if isfield(handles,'hMain')
    Data_Main = guidata(handles.hMain);
    Data_Main.Structure = Structure;
    guidata(handles.hMain,Data_Main)
    
    % change Name of Main GUI to help identifying which Structural Model is
    % using
    Model_Name    = handles.hModel.Name;
    handles.hMain.Name = ['COSMOSS: ' Model_Name];
end

handles.Structure = Structure;
guidata(hObject,handles)

function PlotMolecule(hObject, eventdata, handles)
PlotXYZfiles_Acid(handles.Structure);

function PlotModes(hObject, eventdata, handles)
Plot_Modes(handles.hModel);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hModel_TCO', handles)
disp('Updated handles exported!')
