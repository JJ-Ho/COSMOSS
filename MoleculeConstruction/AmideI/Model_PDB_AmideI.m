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
gui_Singleton = 1;
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
StrucGUI = GUI_PDB_AmideI(hObject);
handles.StrucGUI = StrucGUI; % export GUI handles to handles

% Get Main function's handles
if nargin > 3
    hMain = varargin{1};
    handles.hMain = hMain;
end
% Update handles structure
guidata(hObject, handles);

% Reset Non-Label Frequency, anharmonicity, and F_min/F_Max to fit amideI mode
Data_Main = guidata(handles.hMain);
hMainGUI  = Data_Main.GUI_Main;

set(hMainGUI.NLFreq ,'String','1644')
set(hMainGUI.LFreq  ,'String','1604')
set(hMainGUI.Anharm ,'String','12')
set(hMainGUI.Beta_NN,'String','0.8')
set(hMainGUI.X_Min  ,'String','1550')
set(hMainGUI.X_Max  ,'String','1700')

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
PDB_Path = [PWD, '/PDB_files/'];

[FilesName,PathName,FilterIndex] = uigetfile({'*.pdb','PDB file'; ...
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


function UpdateStructure(hObject, eventdata, handles)

StrucGUI  = handles.StrucGUI;
Data_Main = guidata(handles.hMain);
hMainGUI  = Data_Main.GUI_Main;

%% Read GUI variables
Phi_D        = str2double(get(StrucGUI.Phi  ,'String'));
Psi_D        = str2double(get(StrucGUI.Psi  ,'String'));
Theta_D      = str2double(get(StrucGUI.Theta,'String'));

NLFreq       = str2double(get(hMainGUI.NLFreq  ,'String'));
Anharm       = str2double(get(hMainGUI.Anharm  ,'String'));

%% Construct molecule
Num_Atoms = handles.Num_Atoms;
XYZ       = handles.XYZ;
AtomName  = handles.AtomName;
FilesName = handles.FilesName;

Structure = GetAmideI(Num_Atoms,XYZ,AtomName,FilesName,...
                      'Phi_D',Phi_D,...
                      'Psi_D',Psi_D,...
                      'Theta_D',Theta_D,...
                      'NLFreq',NLFreq,...
                      'Anharm',Anharm);

%% Export result to Main guidata
Data_Main.Structure = Structure;
guidata(handles.hMain,Data_Main)


function PlotMolecule(hObject, eventdata, handles)
Data_Main = guidata(handles.hMain);
PlotXYZfiles_AmideI(Data_Main.Structure)





