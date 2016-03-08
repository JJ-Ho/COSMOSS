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

% Last Modified by GUIDE v2.5 20-Feb-2016 16:34:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
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
GUI_Struc = GUI_Betasheet_AmideI(hObject);
handles.GUI_Struc = GUI_Struc; % export GUI handles to handles

% Get Main function's handles
% Reset Non-Label Frequency, anharmonicity, and F_min/F_Max to fit amideI mode
% check if run this GUI stand along

% Get Main function's handles
if nargin > 3    
    if ishandle(varargin{1}) 
       hMain = varargin{1};
       Data_Main = guidata(hMain);

       handles.hMain = hMain;
       handles.Data_Main = Data_Main;
       
       % PRE ASSIGN VALUES TO SUBSTITUTE MAIN GUI VALUES
       GUI_Main  = Data_Main.GUI_Main;
       set(GUI_Main.Beta_NN,'String','0.8')
       set(GUI_Main.X_Min  ,'String','1550')
       set(GUI_Main.X_Max  ,'String','1700')
    end
else
    disp('Running Model_Betasheet_AmideI in stand alone mode.')    
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Model_Betasheet_AmideI wait for user response (see UIRESUME)
% uiwait(handles.hModel);


% --- Outputs from this function are returned to the command line.
function varargout = Model_Betasheet_AmideI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function UpdateStructure(hObject, eventdata, handles)
%% Read GUI variables
GUI_Struc  = handles.GUI_Struc;
GUI_Inputs = ParseGUI_Betasheet(GUI_Struc);

%% Construct molecule
BB        = ConstuctBetaSheet(GUI_Inputs);
Structure = GetAmideI(BB.Num_Atoms,...
                      BB.XYZ,...
                      BB.AtomName,...
                      BB.FilesName,...
                      GUI_Inputs);

%% Export extra info into Structure                   
Structure.N_Residue         = BB.N_Residue;
Structure.N_Strand          = BB.N_Strand;
Structure.N_Mode_per_Starnd = BB.N_Residue-1;

% C terminus Index
Structure.Ind_H = BB.Ind_H;
Structure.Ind_O = BB.Ind_O;

% Betasheet orientation info export
Structure.TransV = BB.TransV;
Structure.TwistV = BB.TwistV;
Structure.RotV   = [GUI_Inputs.Phi_D,GUI_Inputs.Psi_D,GUI_Inputs.Theta_D];

% Export into Structure so it can be passsed around different GUIs
Structure.StructModel = 4;

%% Export result to Main guidata
% check if this program run stand along
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

disp('Structure file generated!')

function hF = PlotMolecule(hObject, eventdata, handles)
% Read GUI variables
GUI_Struc  = handles.GUI_Struc;
GUI_Inputs = ParseGUI_Betasheet(GUI_Struc);

% Read the Molecule frame to Lab frame orientation from COSMOSS
hMain = handles.hMain;
GUI_Data_Main = guidata(hMain);
GUI_Inputs_Main = ParseGUI_Main(GUI_Data_Main);
% Pass the MF-LB Eular angles to Plotting function
GUI_Inputs.Avg_Phi   = GUI_Inputs_Main.Avg_Phi;
GUI_Inputs.Avg_Theta = GUI_Inputs_Main.Avg_Theta;
GUI_Inputs.Avg_Psi   = GUI_Inputs_Main.Avg_Psi;

hF = Plot_Betasheet_AmideI(handles.Structure,GUI_Inputs);

function PlotModes(hObject, eventdata, handles)
Plot_Modes(handles.hModel);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hModel_Betasheet_AmideI', handles)
disp('Updated handles exported!')
