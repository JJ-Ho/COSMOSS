function varargout = Model_Betasheet_AmideI(varargin)
% MODEL_BETASHEET_AMIDEI MATLAB code for Model_Betasheet_AmideI.fig
%      MODEL_BETASHEET_AMIDEI, by itself, creates a new MODEL_BETASHEET_AMIDEI or raises the existing
%      singleton*.
%
%      H = MODEL_BETASHEET_AMIDEI returns the handle to a new MODEL_BETASHEET_AMIDEI or the handle to
%      the existing singleton*.
%
%      MODEL_BETASHEET_AMIDEI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_BETASHEET_AMIDEI.M with the given input arguments.
%
%      MODEL_BETASHEET_AMIDEI('Property','Value',...) creates a new MODEL_BETASHEET_AMIDEI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_Betasheet_AmideI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_Betasheet_AmideI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Model_Betasheet_AmideI

% Last Modified by GUIDE v2.5 16-Sep-2015 10:37:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Model_Betasheet_AmideI_OpeningFcn, ...
                   'gui_OutputFcn',  @Model_Betasheet_AmideI_OutputFcn, ...
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


% --- Executes just before Model_Betasheet_AmideI is made visible.
function Model_Betasheet_AmideI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_Betasheet_AmideI (see VARARGIN)

% Choose default command line output for Model_Betasheet_AmideI
handles.output = hObject;

% Call createInterface to create GUI elements
StrucGUI = GUI_Betasheet_AmideI(hObject);
handles.StrucGUI = StrucGUI; % export GUI handles to handles

% Get Main function's handles
if nargin > 3
    hMain = varargin{1};
    handles.hMain = hMain;
end

% Update handles structure
guidata(hObject, handles);

% Reset Non-Label Frequency, anharmonicity, and F_min/F_Max to fit amideI mode
% check if run this GUI stand along
if isfield(handles,'hMain')
    Data_Main = guidata(handles.hMain);
    hMainGUI  = Data_Main.GUI_Main;

    set(hMainGUI.NLFreq ,'String','1644')
    set(hMainGUI.LFreq  ,'String','1604')
    set(hMainGUI.Anharm ,'String','12')
    set(hMainGUI.Beta_NN,'String','0.8')
    set(hMainGUI.X_Min  ,'String','1550')
    set(hMainGUI.X_Max  ,'String','1700')
else 
    disp('Running in stand alone mode.')
end
% UIWAIT makes Model_Betasheet_AmideI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Model_Betasheet_AmideI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function UpdateStructure(hObject, eventdata, handles)

StrucGUI  = handles.StrucGUI;
%% Read GUI variables
N_Residue = str2double(get(StrucGUI.N_Residue,'String'));
N_Strand  = str2double(get(StrucGUI.N_Strand ,'String'));
Trans_X   = str2double(get(StrucGUI.Trans_X  ,'String'));
Trans_Y   = str2double(get(StrucGUI.Trans_Y  ,'String'));
Trans_Z   = str2double(get(StrucGUI.Trans_Z  ,'String'));
Rot_X     = str2double(get(StrucGUI.Rot_X    ,'String'));
Rot_Y     = str2double(get(StrucGUI.Rot_Y    ,'String'));
Rot_Z     = str2double(get(StrucGUI.Rot_Z    ,'String'));

% check if run this GUI stand along
if isfield(handles,'hMain')
    Data_Main = guidata(handles.hMain);
    hMainGUI  = Data_Main.GUI_Main;
    
    NLFreq    = str2double(get(hMainGUI.NLFreq  ,'String'));
    Anharm    = str2double(get(hMainGUI.Anharm  ,'String'));
else
    NLFreq = 1644;
    Anharm = 12;
end

%% Construct molecule
TransV = [Trans_X,Trans_Y,Trans_Z];
RotV   = [Rot_X,Rot_Y,Rot_Z];

XYZ       = ConstuctBetaSheet(N_Residue,N_Strand,TransV,RotV);
Structure = GetAmideI_CON_XYZ_betasheet(XYZ);

%% Export result to Main guidata
Data_Main.Structure = Structure;
% check if this program run stand along
if isfield(handles,'hMain')
    guidata(handles.hMain,Data_Main)
else
    handles.Structure = Structure;
    guidata(hObject,handles)
end

disp('Structure file generated!')


function PlotMolecule(hObject, eventdata, handles)

% check if this program run stand along
if isfield(handles,'hMain')
    Data_Main = guidata(handles.hMain);
    Plot_Betasheet_AmideI(Data_Main.Structure)
else
    Plot_Betasheet_AmideI(handles.Structure)
end