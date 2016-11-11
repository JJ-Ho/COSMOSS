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

% Last Modified by GUIDE v2.5 03-Nov-2016 18:18:13

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

% Get Main (upper level) function's handles
if nargin > 3    
    if ishandle(varargin{1}) 
       hMain = varargin{1};
       Data_Main = guidata(hMain);

       handles.hMain = hMain;
       handles.Data_Main = Data_Main;
    else
        if strcmp(varargin{1},'Comb2_Mode')
            disp('Running Model_TwoDGrid as a sub GUI of Comb2')
        end
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
% varargout{1} = handles.output;
varargout{1} = handles; % export the whole guidata instead of GUI base figure handle.

function LoadG09(hObject, eventdata, handles)
%% retreive GUI handles
GUI_Struc  = handles.GUI_Struc;

%% Call uigetfile for G09 path
PWD = pwd;
G09_default_folder = [PWD, '/StructureFiles/G09/'];

[FilesName,PathName,~] = uigetfile({'*.txt','Formatted G09 output'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',G09_default_folder);
    
G09_Path = [PathName FilesName];                             
disp([FilesName,' loaded...'])

%% Load selected G09 file
G09_Output = ReadG09Input(G09_Path);
G09_Output.G09_Path = G09_Path;

%% Export to Model handles
handles.Structure.G09_Output = G09_Output;
guidata(hObject,handles)

%% update structure and the file name on GUI
handles.hModel.Name = ['Model_2DGrid: ', FilesName];
% update mol frame Atom index
GUI_Struc.MF_Center.String = num2str(G09_Output.Mol_Frame.Center_Ind);
GUI_Struc.MF_Zi.String     = num2str(G09_Output.Mol_Frame.Z_i_Ind);
GUI_Struc.MF_Zf.String     = num2str(G09_Output.Mol_Frame.Z_f_Ind);
GUI_Struc.MF_XYi.String    = num2str(G09_Output.Mol_Frame.XY_i_Ind);
GUI_Struc.MF_XYf.String    = num2str(G09_Output.Mol_Frame.XY_f_Ind);

UpdateStructure(hObject, eventdata, handles)

function UpdateStructure(hObject, eventdata, handles)
%% retreive GUI inputs
GUI_Struc  = handles.GUI_Struc;
GUI_Inputs = ParseGUI_TwoDGrid(GUI_Struc);
G09_Output = handles.Structure.G09_Output;

% update the mol frame info into G09_Output
GUI_MF_Info.Center_Ind = GUI_Inputs.MF_Center;
GUI_MF_Info.Z_i_Ind    = GUI_Inputs.MF_Zi;
GUI_MF_Info.Z_f_Ind    = GUI_Inputs.MF_Zf;
GUI_MF_Info.XY_i_Ind   = GUI_Inputs.MF_XYi;
GUI_MF_Info.XY_f_Ind   = GUI_Inputs.MF_XYf;
GUI_MF_Info.Frame_Type = GUI_Inputs.Frame_Type;
GUI_MF_Info.BondAvg    = GUI_Inputs.BondAvg;

GUI_LF_Info.LF_Phi   = GUI_Inputs.LF_Phi./180*pi;
GUI_LF_Info.LF_Psi   = GUI_Inputs.LF_Psi./180*pi;
GUI_LF_Info.LF_Theta = GUI_Inputs.LF_Theta./180*pi;

% update the molecule frame info into G09_Output
G09_Output.MF_Info = GUI_MF_Info;
G09_Output.LF_Info = GUI_LF_Info;

%% Rot/Trans the raw ouput of G09 to molecule frame ane deal with rotatable bond average 
% Rotate the structre from G09 simulation frame to defined molecule frame 
% then do bond rotational average if selected. Finally rotate molecule into 
% lab frame so the monomer is ready to form 2D grid.
Monomer = R_GF2LF(G09_Output);

%% Construct 2D grid
Structure = ConstructGrid(Monomer.LF,GUI_Inputs);

% Export into Structure so it can be passsed around different GUIs
Structure.StructModel = 3;

%% Export result to Main guidata
% include G09 ouput to structure
Structure.G09_Output = G09_Output;

% include FieldName of GUI Inputs
[~,~,~,hGUIParser] = StructureModel(Structure.StructModel);
[~,GUI_FieldName] = hGUIParser(GUI_Struc);
handles.GUI_FieldName = GUI_FieldName;

% update handles
handles.Structure  = Structure;
handles.GUI_Inputs = GUI_Inputs;
guidata(hObject,handles)

% update to other GUIs
Export2GUIs(handles)

disp('Structure file generated!')


function hF = PlotMolecule(hObject, eventdata, handles)
%% Update structure
UpdateStructure(hObject, eventdata, handles)

%% Read GUI
GUI_Struc  = handles.GUI_Struc;
GUI_Inputs = ParseGUI_TwoDGrid(GUI_Struc);

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

hF = PlotXYZ_Grid(handles.Structure,GUI_Inputs);

function PlotModes(hObject, eventdata, handles)
Plot_Modes(handles.hModel);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hModel_TwoDGrid', handles)
disp('Updated handles exported!')

function Replace_Convention_Label(hObject, eventdata, handles)
% retreive GUI inputs
GUI_Struc  = handles.GUI_Struc;
Frame_Type = GUI_Struc.Frame_Type.Value;

switch Frame_Type
    case 1 %Frame_Type = 'XZ';
        GUI_Struc.MF_XYi_text.String = 'XZ_i Ind.:';
        GUI_Struc.MF_XYf_text.String = 'XZ_f Ind.:';
    case 2 %Frame_Type = 'YZ';
        GUI_Struc.MF_XYi_text.String = 'YZ_i Ind.:';
        GUI_Struc.MF_XYf_text.String = 'YZ_f Ind.:';
end