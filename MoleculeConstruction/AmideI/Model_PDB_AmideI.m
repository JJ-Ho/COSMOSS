% ^ GUI Skeleton ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

function Model_PDB_AmideI_OpeningFcn(hModel_PDB_AmideI, eventdata, GUI_data, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_PDB_AmideI (see VARARGIN)
if nargin > 3
    switch varargin{1}
        case 'COSMOSS'
            hCOSMOSS = varargin{2};
            Data_COSMOSS = guidata(hCOSMOSS);
            
            %PRE ASSIGN VALUES TO SUBSTITUTE MAIN GUI VALUES
            hGUIs_COSMOSS  = Data_COSMOSS.hGUIs;
            set(hGUIs_COSMOSS.Beta_NN,'String','0.8')
            set(hGUIs_COSMOSS.X_Min  ,'String','1550')
            set(hGUIs_COSMOSS.X_Max  ,'String','1750')
            
            GUI_data.hCOSMOSS = hCOSMOSS;
            
            disp('Running Model_PDB_AmideI durectly from COSMOSS...')
        case 'Comb2'
            hModel_Comb2 = varargin{2};
            Comb2_Order  = varargin{3};
            
            GUI_data.hModel_Comb2 = hModel_Comb2;
            GUI_data.Comb2_Order  = Comb2_Order;
            
            % Add comb2 order # to GUI title, if necessary
            TitleStr = hModel_PDB_AmideI.Name;
            if ~strcmp(TitleStr(1),'#')
                hModel_PDB_AmideI.Name = ['#',int2str(Comb2_Order),', ',TitleStr];
            end
            
            disp('Running Model_PDB_AmideI as a sub GUI of Comb2...')
    end
else
    disp('Running Model_PDB_AmideI in stand alone mode...')
end

% Call create Interface to create GUI elements
hGUIs = GUI_PDB_AmideI(hModel_PDB_AmideI);

% Prep necessary data to be export
GUI_data.hModel_PDB_AmideI = hModel_PDB_AmideI;
GUI_data.hGUIs             = hGUIs;

% Update handles structure
guidata(hModel_PDB_AmideI, GUI_data);

% --- Outputs from this function are returned to the command line.
function varargout = Model_PDB_AmideI_OutputFcn(hModel_PDB_AmideI, GUI_data, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hModel_PDB_AmideI;

function Export_Handle_Callback(hObject, eventdata, GUI_data)
% export handles back to work space
assignin('base', 'Data_PDB_AmideI', GUI_data)
disp('Updated GUI Data_PDB_AmideI exported!')
% ^ GUI Skeleton ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    



function LoadStructure(hObject, eventdata, GUI_data)
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
GUI_data.Num_Atoms = Num_Atoms;
GUI_data.XYZ       = XYZ;
GUI_data.AtomName  = AtomName;
GUI_data.FilesName = FilesName;

guidata(hObject,GUI_data)

%% Update PDB name and data to GUI
set(GUI_data.hGUIs.PDB_Name,'String',FilesName)
UpdateStructure(hObject, eventdata, GUI_data)

function UpdateStructure(hObject, eventdata, GUI_data)
%% Read GUI variables
hGUIs  = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_AmideI(hGUIs);
    
%% Construct molecule
Num_Atoms = GUI_data.Num_Atoms;
XYZ       = GUI_data.XYZ;
AtomName  = GUI_data.AtomName;
FilesName = GUI_data.FilesName;

Structure = GetAmideI(Num_Atoms,XYZ,AtomName,FilesName,GUI_Inputs);
                  
% Export into Structure so it can be passsed around different GUIs
Structure.StructModel = 2;

%% Export result to Main guidata
GUI_data.Structure = Structure;

% include FieldName of GUI Inputs
[~,~,~,hGUIParser] = StructureModel(Structure.StructModel);
[~,GUI_FieldName] = hGUIParser(hGUIs);
GUI_data.GUI_FieldName = GUI_FieldName;

guidata(hObject,GUI_data)

% update to other GUIs
Export2GUIs(GUI_data)

disp('Structure file generated!')

function hF = PlotMolecule(hObject, eventdata, GUI_data)
% Read GUI variables
hGUIs  = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_AmideI(hGUIs);

%- This part is obsolete, since the lab frame ensemble avg should not take
%  orientation inputs, will be removed later
% Read the Molecule frame to Lab frame orientation from COSMOSS
% hMain = handles.hMain;
% GUI_Data_Main = guidata(hMain);
% GUI_Inputs_Main = ParseGUI_Main(GUI_Data_Main);
% % Pass the MF-LB Eular angles to Plotting function
% GUI_Inputs.Avg_Phi   = GUI_Inputs_Main.Avg_Phi;
% GUI_Inputs.Avg_Theta = GUI_Inputs_Main.Avg_Theta;
% GUI_Inputs.Avg_Psi   = GUI_Inputs_Main.Avg_Psi;

GUI_Inputs.Avg_Phi   = 0;
GUI_Inputs.Avg_Theta = 0;
GUI_Inputs.Avg_Psi   = 0;
%--------------------------------------------------------------------------


hF = PlotXYZfiles_AmideI(GUI_data.Structure,GUI_Inputs);

function PlotModes(hObject, eventdata, GUI_data)
Plot_Modes(GUI_data.hModel_PDB_AmideI);






