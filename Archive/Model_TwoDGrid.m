%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

function Model_TwoDGrid_OpeningFcn(hModel_TwoDGrid, eventdata, GUI_data, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_TwoDGrid (see VARARGIN)
if nargin > 3
    switch varargin{1}
        case 'COSMOSS'
            hCOSMOSS = varargin{2};
            Data_COSMOSS = guidata(hCOSMOSS);
            
            %PRE ASSIGN VALUES TO SUBSTITUTE MAIN GUI VALUES
            hGUIs_COSMOSS  = Data_COSMOSS.hGUIs;
            set(hGUIs_COSMOSS.F_Min  ,'String','1650')
            set(hGUIs_COSMOSS.F_Max  ,'String','1750')
            
            GUI_data.hCOSMOSS = hCOSMOSS;
            
            disp('Running Model_TwoDGrid directly from COSMOSS...')
        case 'Comb2'
            hModel_Comb2 = varargin{2};
            Comb2_Order  = varargin{3};
            
            GUI_data.hModel_Comb2 = hModel_Comb2;
            GUI_data.Comb2_Order  = Comb2_Order;
            
            % Add comb2 order # to GUI title, if necessary
            TitleStr = hModel_TwoDGrid.Name;
            if ~strcmp(TitleStr(1),'#')
                hModel_TwoDGrid.Name = ['#',int2str(Comb2_Order),', ',TitleStr];
            end
            
            disp('Running Model_TwoDGrid as a sub GUI of Comb2...')
    end
else
    disp('Running Model_TwoDGrid in stand alone mode...')
end

% Call create Interface to create GUI elements
hGUIs = GUI_TwoDGrid(hModel_TwoDGrid);

% Prep necessary data to be saved in GUI_data
GUI_data.hModel_TwoDGrid = hModel_TwoDGrid;
GUI_data.hGUIs           = hGUIs; % export GUI handles to handles

guidata(hModel_TwoDGrid,GUI_data);

function varargout = Model_TwoDGrid_OutputFcn(hModel_TwoDGrid, eventdata, GUI_data) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = hModel_TwoDGrid; % export the whole guidata instead of GUI base figure handle.

function Export_Handle_Callback(hObject, eventdata, GUI_data)
% export handles back to work space
assignin('base', 'Data_TwoDGrid', GUI_data)
disp('Updated GUI Data_TwoDGrid exported!')
%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



function LoadG09(hObject, eventdata, GUI_data)
%% retreive GUI handles
hGUIs  = GUI_data.hGUIs;

%% Call uigetfile for G09 path
PWD = pwd;
G09_default_folder = [PWD, '/StructureFiles/G09/'];

[FilesName,PathName,~] = uigetfile({'*.txt','Formatted G09 output'; ...
                                    '*,*','All Files'},...
                                    'Select inputs',G09_default_folder);
    
G09_Path = [PathName FilesName];                             
disp([FilesName,' loaded...'])

%% Load selected G09 file
S_G09_MF = ReadG09Input(G09_Path);
% G09_Output.G09_Path = G09_Path;

% Add function handles for plotting
S_G09_MF.hPlotFunc = @PlotXYZ_Grid;
S_G09_MF.hParseGUIFunc = @ParseGUI_TwoDGrid;
S_G09_MF.hGUIs = hGUIs;

%% Export to Model handles
GUI_data.G09 = S_G09_MF;
guidata(hObject,GUI_data)

%% update structure and the file name on GUI
TitleStr = '2DGrid: '; 
if isfield(GUI_data,'Comb2_Order')
    TitleStr = ['#',int2str(GUI_data.Comb2_Order),', ',TitleStr];
end
GUI_data.hModel_TwoDGrid.Name = [TitleStr, FilesName];

% update G09 info onto GUI
hGUIs.MF_Center.String = num2str(S_G09_MF.Extra.Mol_Frame.Center_Ind);
hGUIs.MF_Zi.String     = num2str(S_G09_MF.Extra.Mol_Frame.Z_i_Ind);
hGUIs.MF_Zf.String     = num2str(S_G09_MF.Extra.Mol_Frame.Z_f_Ind);
hGUIs.MF_XZi.String    = num2str(S_G09_MF.Extra.Mol_Frame.XZ_i_Ind);
hGUIs.MF_XZf.String    = num2str(S_G09_MF.Extra.Mol_Frame.XZ_f_Ind);
hGUIs.NLFreq.String    = num2str(S_G09_MF.LocFreq);
hGUIs.Anharm.String    = num2str(S_G09_MF.LocAnharm);

UpdateStructure(hObject, eventdata, GUI_data)

function UpdateStructure(hObject, eventdata, GUI_data)
%% retreive GUI inputs
hGUIs  = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_TwoDGrid(hGUIs);
G09 = GUI_data.G09;

%% Rot/Trans the raw ouput of G09 to molecule frame ane deal with rotatable bond average 
% Rotate the structre from G09 simulation frame to defined molecule frame 
% then do bond rotational average if selected. Finally rotate molecule into 
% lab frame so the monomer is ready to form 2D grid.
S_Monomer_LF = R_MF2LF(G09,GUI_Inputs);

%% Construct 2D grid
Structure_Grid = ConstructGrid(S_Monomer_LF,GUI_Inputs);

% Export into Structure so it can be passsed around different GUIs
Structure_Grid.StructModel = 3;

% export necessary handle and functions
Structure_Grid.hPlotFunc = @PlotXYZ_Grid;
Structure_Grid.hParseGUIFunc = @ParseGUI_TwoDGrid;
Structure_Grid.hGUIs = hGUIs;

%% Export result to Main guidata
% include FieldName of GUI Inputs
[~,~,~,fhGUIParser] = StructureModel(Structure_Grid.StructModel);
[~,GUI_FieldName] = fhGUIParser(hGUIs);
GUI_data.GUI_FieldName = GUI_FieldName;

% update handles
GUI_data.Monomer    = S_Monomer_LF;
GUI_data.Structure  = Structure_Grid;
GUI_data.GUI_Inputs = GUI_Inputs;
guidata(hObject,GUI_data)

% update to other GUIs
Export2GUIs(GUI_data)

disp('Structure file generated!')

function hF = PlotMolecule(hObject, eventdata, GUI_data)
hF = GUI_data.Structure.Draw;

function PlotModes(hObject, eventdata, GUI_data)
Plot_Modes(GUI_data.hModel_TwoDGrid);

function Replace_Convention_Label(hObject, eventdata, GUI_data)
% retreive GUI inputs
hGUIs  = GUI_data.hGUIs;
Frame_Type = hGUIs.Frame_Type.Value;

switch Frame_Type
    case 1 %Frame_Type = 'XZ';
        hGUIs.MF_XYi_text.String = 'XZ_i Ind.:';
        hGUIs.MF_XYf_text.String = 'XZ_f Ind.:';
    case 2 %Frame_Type = 'YZ';
        hGUIs.MF_XYi_text.String = 'YZ_i Ind.:';
        hGUIs.MF_XYf_text.String = 'YZ_f Ind.:';
end