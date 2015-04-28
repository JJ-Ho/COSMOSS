function varargout = Model_TwoDGrid(varargin)
% MODEL_TWODGRID MATLAB code for Model_TwoDGrid.fig
%      MODEL_TWODGRID, by itself, creates a new MODEL_TWODGRID or raises the existing
%      singleton*.
%
%      H = MODEL_TWODGRID returns the handle to a new MODEL_TWODGRID or the handle to
%      the existing singleton*.
%
%      MODEL_TWODGRID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_TWODGRID.M with the given input arguments.
%
%      MODEL_TWODGRID('Property','Value',...) creates a new MODEL_TWODGRID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_TwoDGrid_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_TwoDGrid_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Model_TwoDGrid

% Last Modified by GUIDE v2.5 27-Apr-2015 21:32:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Model_TwoDGrid_OpeningFcn, ...
                   'gui_OutputFcn',  @Model_TwoDGrid_OutputFcn, ...
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


% --- Executes just before Model_TwoDGrid is made visible.
function Model_TwoDGrid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_TwoDGrid (see VARARGIN)

% Choose default command line output for Model_TwoDGrid
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Struc = GUI_TwoDGrid(hObject);
handles.GUI_Struc = GUI_Struc; % export GUI handles to handles

% Get Main function's handles
if nargin > 3
    hMain = varargin{1};
    Data_Main = guidata(handles.hMain);
else
    hMain = 'Debug';
    Data_Main.GUI_Main.NLFreq = 1720;
end
handles.hMain = hMain;
handles.Data_Main = Data_Main;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Model_TwoDGrid wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Model_TwoDGrid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function LoadG09(hObject, eventdata, handles)

function UpdateStructure(hObject, eventdata, handles)
GUI_Struc = handles.GUI_Struc;
Data_Main = handles.Data_Main;
GUI_Main  = Data_Main.GUI_Main;

%% Construct molecule
GUI_Inputs = ParseGUI_TwoDGrid(GUI_Struc);
% GUI_Inputs.NLFreq = str2double(get(GUI_Main.NLFreq  ,'String'));

Structure  = ConstructGrid(GUI_Inputs);

%% Export result to Main guidata
Data_Main.Structure = Structure;
% guidata(handles.hMain,Data_Main)


function PlotMolecule(hObject, eventdata, handles)
