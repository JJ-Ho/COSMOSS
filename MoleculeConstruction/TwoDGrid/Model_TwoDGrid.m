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

% Last Modified by GUIDE v2.5 10-Mar-2016 15:45:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
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
    if ishandle(varargin{1}) 
       hMain = varargin{1};
       Data_Main = guidata(hMain);

       handles.hMain = hMain;
       handles.Data_Main = Data_Main;
    end
else
    disp('Running Model_TwoDGrid in stand alone mode.')    
end

% Update handles structure
guidata(hObject,handles);

% UIWAIT makes Model_TwoDGrid wait for user response (see UIRESUME)
% uiwait(handles.hModel);


% --- Outputs from this function are returned to the command line.
function varargout = Model_TwoDGrid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function LoadG09(hObject, eventdata, handles)
%% Call uigetfile for G09 path
PWD = pwd;
G09_default_folder = [PWD, '/StructureFiles/G09/'];

[FilesName,PathName,~] = uigetfile({'*.txt','Formatted G09 output'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',G09_default_folder);
    
G09_Path = [PathName FilesName];                             
disp([FilesName,' loaded...'])

%% Export to Model handles
handles.G09_Path = G09_Path;
guidata(hObject,handles)

%% update structure and the file name on GUI
set(handles.GUI_Struc.G09_FileName,'String',['G09 file: ', FilesName])
UpdateStructure(hObject, eventdata, handles)

function UpdateStructure(hObject, eventdata, handles)
% retreive GUI inputs
GUI_Struc  = handles.GUI_Struc;
GUI_Inputs = ParseGUI_TwoDGrid(GUI_Struc);

%% Load selected G09 file
Ang_Phi      = GUI_Inputs.Ang_Phi;
Ang_Psi      = GUI_Inputs.Ang_Psi;
Ang_Theta    = GUI_Inputs.Ang_Theta;
Eular_MF_D   = [Ang_Phi,Ang_Psi,Ang_Theta];
Eular_MF_R   = Eular_MF_D./180*pi; % turn to radius unit

Monomer_Axes = GUI_Inputs.Monomer_Axes;
switch Monomer_Axes
    case 1
        MolFrame_Convention = 'XZ';
    case 2
        MolFrame_Convention = 'YZ';
end
Gaussian_Input = ReadG09Input(handles.G09_Path,'MolFrame',MolFrame_Convention);

%% Construct molecule
Structure = ConstructGrid(Gaussian_Input,GUI_Inputs);

% Export into Structure so it can be passsed around different GUIs
Structure.StructModel = 3;

%% Export result to Main guidata

% check if this program run stand along, if not push Structure info to Main
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

function hF = PlotMolecule(hObject, eventdata, handles)
%% Update structure
UpdateStructure(hObject, eventdata, handles)

%% Read GUI
GUI_Struc  = handles.GUI_Struc;
GUI_Inputs = ParseGUI_TwoDGrid(GUI_Struc);

% Read the Molecule frame to Lab frame orientation from COSMOSS
hMain = handles.hMain;
GUI_Data_Main = guidata(hMain);
GUI_Inputs_Main = ParseGUI_Main(GUI_Data_Main);
% Pass the MF-LB Eular angles to Plotting function
GUI_Inputs.Avg_Phi   = GUI_Inputs_Main.Avg_Phi;
GUI_Inputs.Avg_Theta = GUI_Inputs_Main.Avg_Theta;
GUI_Inputs.Avg_Psi   = GUI_Inputs_Main.Avg_Psi;

hF = PlotXYZ_Grid(handles.Structure,GUI_Inputs);

function PlotModes(hObject, eventdata, handles)
Plot_Modes(handles.hModel);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hModel_TwoDGrid', handles)
disp('Updated handles exported!')
