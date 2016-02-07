function varargout = Model_PDB_AmideI(varargin)
% Model_PDB_AmideI MATLAB code for Model_PDB_AmideI.fig
%      Model_PDB_AmideI, by itself, creates a new Model_PDB_AmideI or raises the existing
%      singleton*.
%
%      H = Model_PDB_AmideI returns the handle to a new Model_PDB_AmideI or the handle to
%      the existing singleton*.
%
%      Model_PDB_AmideI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Model_PDB_AmideI.M with the given input arguments.
%
%      Model_PDB_AmideI('Property','Value',...) creates a new Model_PDB_AmideI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_PDB_AmideI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_PDB_AmideI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Model_PDB_AmideI

% Last Modified by GUIDE v2.5 01-Oct-2014 16:16:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Model_PDB_AmideI_OpeningFcn, ...
                   'gui_OutputFcn',  @Model_PDB_AmideI_OutputFcn, ...
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


% --- Executes just before Model_PDB_AmideI is made visible.
function Model_PDB_AmideI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_PDB_AmideI (see VARARGIN)

% Choose default command line output for Model_PDB_AmideI
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Struc = GUI_PDB_AmideI(hObject);
handles.GUI_Struc = GUI_Struc; % export GUI handles to handles

% Reset Non-Label Frequency, anharmonicity, and F_min/F_Max to fit amideI mode
% check if run this GUI stand along
if nargin > 3    
    if ishandle(varargin{1}) 
        hMain = varargin{1};
        Data_Main = guidata(hMain);
        
        handles.hMain = hMain;
        handles.Data_Main = Data_Main;
        
        % PRE ASSIGN VALUES TO SUBSTITUTE MAIN GUI VALUES
        GUI_Main  = Data_Main.GUI_Main;
        set(GUI_Main.NLFreq ,'String','1644')
        set(GUI_Main.LFreq  ,'String','1604')
        set(GUI_Main.Anharm ,'String','12')
        set(GUI_Main.Beta_NN,'String','0.8')
        set(GUI_Main.X_Min  ,'String','1550')
        set(GUI_Main.X_Max  ,'String','1700')
    end
else
    disp('Running Model_PDB_AmideI in stand alone mode.')    
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Model_PDB_AmideI wait for user response (see UIRESUME)
% uiwait(handles.Model_PDB_AmideI);


% --- Outputs from this function are returned to the command line.
function varargout = Model_PDB_AmideI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function LoadStructure(hObject, eventdata, handles)
%% Get pdb file location 
PWD = pwd;
PDB_Path = [PWD, '/StructureFiles/PDB/'];

[FilesName,PathName,~] = uigetfile({'*.pdb','PDB file'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',PDB_Path);

%% Parse molecule structure

PDB = pdbread([PathName FilesName]);
Atom_Data = PDB.Model.Atom;
Num_Atoms = size(Atom_Data,2);

% Get coordination data
XYZ = zeros(Num_Atoms,3);
AtomName = cell(Num_Atoms,1);
for II = 1:Num_Atoms
    A = Atom_Data(II);
    XYZ(II,:) = [A.X, A.Y, A.Z];
    AtomName{II} = Atom_Data(II).AtomName;
end

%% output to GUI
handles.Num_Atoms = Num_Atoms;
handles.XYZ       = XYZ;
handles.AtomName  = AtomName;
handles.FilesName = FilesName;

guidata(hObject,handles)

%% Update PDB name and data to GUI
set(handles.GUI_Struc.PDB_Name,'String',FilesName)
UpdateStructure(hObject, eventdata, handles)



function UpdateStructure(hObject, eventdata, handles)
%% Read GUI variables
GUI_Struc  = handles.GUI_Struc;
GUI_Inputs = ParseGUI_AmideI(GUI_Struc);
    
%% Construct molecule
Num_Atoms = handles.Num_Atoms;
XYZ       = handles.XYZ;
AtomName  = handles.AtomName;
FilesName = handles.FilesName;

Structure = GetAmideI(Num_Atoms,XYZ,AtomName,FilesName,GUI_Inputs);
                  
% Export into Structure so it can be passsed around different GUIs
Structure.StructModel = 2;

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

function PlotMolecule(hObject, eventdata, handles)
PlotXYZfiles_AmideI(handles.Structure)

function PlotModes(hObject, eventdata, handles)
Plot_Modes(handles.hModel);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hModel_PDB_AmideI', handles)
disp('Updated handles exported!')
    




