function varargout = Plot_Modes(varargin)
% PLOT_MODES MATLAB code for Plot_Modes.fig
%      PLOT_MODES, by itself, creates a new PLOT_MODES or raises the existing
%      singleton*.
%
%      H = PLOT_MODES returns the handle to a new PLOT_MODES or the handle to
%      the existing singleton*.
%
%      PLOT_MODES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_MODES.M with the given input arguments.
%
%      PLOT_MODES('Property','Value',...) creates a new PLOT_MODES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Plot_Modes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Plot_Modes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Plot_Modes

% Last Modified by GUIDE v2.5 10-Feb-2016 16:58:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Plot_Modes_OpeningFcn, ...
                   'gui_OutputFcn',  @Plot_Modes_OutputFcn, ...
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

% --- Executes just before Plot_Modes is made visible.
function Plot_Modes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Plot_Modes (see VARARGIN)

% Choose default command line output for Plot_Modes
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Modes = GUI_Plot_Modes(hObject);

% Get Structural modeling GUI's handles
if nargin > 3    
    if ishandle(varargin{1}) 
       hModel = varargin{1};
    end
else
    disp('Running in stand alone mode, using TCO modle for debug purpose')  
    hModel = Model_TCO(handles);
end

% Change Names on Plot_Exciton GUI to identify which Structural model is
% using
Model_Name = hModel.Name;
Plot_Modes_GUI_Name = GUI_Modes.hPlot_Modes.Name;
GUI_Modes.hPlot_Modes.Name = [Plot_Modes_GUI_Name, ': ', Model_Name];

% Export the handle of Structure Modle to guidata of Plot_Exciton
handles.hModel = hModel;
handles.GUI_Modes = GUI_Modes; % export GUI handles to handles
guidata(hObject, handles);

% update exciton info in Plot_Exciton
Update_Modes(hObject, eventdata, handles)

% --- Outputs from this function are returned to the command line.
function varargout = Plot_Modes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Update_Modes(hObject, eventdata, handles)
% Run OneDSFG to get the corresponding mu and alpha of exciton modes
% Retrieve Label index and Coupling model from COSMOSS GUI if any, if
% running Plot_Exciton stand alone for debug testing, give a field 'debug'
% to use the default values in OneDSFG_Main.m
GUI_Data_hModel = guidata(handles.hModel);
if isfield(GUI_Data_hModel,'hMain')
    GUI_Data_hMain = guidata(GUI_Data_hModel.hMain);
    MainGUI_Inputs = ParseGUI_Main(GUI_Data_hMain);
    %disp('Plot_Modes: Using coupling info from Main GUI')
else
    MainGUI_Inputs.debug = 'debug';
    GUI_Data_hModel.hMain = 'debug';
    %disp('Plot_Modes: Coupling info comes from defulat setting of OneDSFG_Main.m')
end

Structure = GUI_Data_hModel.Structure;
Modes     = Update_Modes_Table(Structure, MainGUI_Inputs);

%% Update handles structure
handles.OneDSFG        = Modes.OneDSFG;
handles.ModeList       = Modes.ModeList;
handles.Structure      = Structure;
handles.hMain          = GUI_Data_hModel.hMain;
% handles.MianGUI_Inputs = MainGUI_Inputs;
guidata(hObject, handles);

% update the list on hPlot_Exciton GUI
set(handles.GUI_Modes.ModeList,'Data',Modes.ModeList)

% call sorting to sort table with the same GUI setting
uitable_SortCallback(hObject, eventdata, handles)

function Update_Figure(hObject, eventdata, handles)
%% Use Update Modes to update the structure and the corresponding Mu & Alpha 
Update_Modes(hObject, eventdata, handles)
handles = guidata(hObject);

GUI_Inputs = ParseGUI_Modes(handles);
Structure  = handles.Structure;
OneDSFG    = handles.OneDSFG;
hModel     = handles.hModel;

% Read the Molecule frame to Lab frame orientation from COSMOSS
hMain = handles.hMain;
GUI_Data_Main = guidata(hMain);
GUI_Inputs_Main = ParseGUI_Main(GUI_Data_Main);
% Pass the MF-LB Eular angles to Plotting function
GUI_Inputs.Avg_Phi   = GUI_Inputs_Main.Avg_Phi;
GUI_Inputs.Avg_Theta = GUI_Inputs_Main.Avg_Theta;
GUI_Inputs.Avg_Psi   = GUI_Inputs_Main.Avg_Psi;

%% Draw molecule by calling the PlotMolecule function in each model
[hFunc_Model,~,~] = StructureModel(Structure.StructModel);
hF = hFunc_Model('PlotMolecule',hModel,eventdata,guidata(hModel));

%% Call Update Figure function
Fig_Output = Update_Modes_Figure(hF, GUI_Inputs, Structure, OneDSFG);

%% update handles
handles.hF      = hF;
handles.Loc_Ind = Fig_Output.Loc_Ind;
handles.Ex_Ind  = Fig_Output.Ex_Ind;
guidata(hObject,handles)

function uitable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
% handles    structure with handles and user data (see GUIDATA)

TableData = handles.GUI_Modes.ModeList.Data;

CurrentCell = eventdata.Indices;
CurrentRowInd = CurrentCell(:,1)';
Mode_Ind_Str = num2str(TableData(CurrentRowInd,1)');

% Update the Mode index on GUI
set(handles.GUI_Modes.Mode_Ind,'String', Mode_Ind_Str);

%% update handles
handles.Mode_Ind_Str = Mode_Ind_Str;
guidata(hObject,handles)

function uitable_SortCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
% handles    structure with handles and user data (see GUIDATA)

TableData = handles.GUI_Modes.ModeList.Data;
SortColumn = handles.GUI_Modes.SortInd.Value;

[~,SortInd] = sort(abs(TableData(:,SortColumn)),'descend');
SortedData = TableData(SortInd,:);

% Update Table on GUI
set(handles.GUI_Modes.ModeList,'Data', SortedData);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hPlot_Modes', handles)
disp('Updated handles exported!')
