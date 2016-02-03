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

% Last Modified by GUIDE v2.5 25-Jan-2016 12:03:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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
handles.GUI_Modes = GUI_Modes; % export GUI handles to handles

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
guidata(hObject, handles);

% update exciton info in Plot_Exciton
Update_Modes(hObject, eventdata, handles)

% UIWAIT makes Plot_Modes wait for user response (see UIRESUME)
% uiwait(handles.hPlot_Modes);


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
    %disp('Using the labeing index and coupling info from Main GUI')
else
    MainGUI_Inputs.debug = 'debug';
    GUI_Data_hModel.hMain = 'debug';
    %disp('Labeling index and coupling info come from defulat setting of OneDSFG_Main.m')
end

Structure = GUI_Data_hModel.Structure;
OneDSFG = OneDSFG_Main(Structure,MainGUI_Inputs);

Ex_Freq     = OneDSFG.H.Sort_Ex_Freq(2:end);
Num_Ex_Mode = length(Ex_Freq);
Ex_Ind      = (1:Num_Ex_Mode)';
Ex_Mu       = squeeze(OneDSFG.Mu.Trans_Ex(1,2:end,:));
Ex_Mu_Z     = Ex_Mu(:,3);
Ex_Mu_Int   = sqrt(sum(Ex_Mu.^2,2));

Ex_Alpha    = squeeze(OneDSFG.Alpha.Trans_Ex(1,2:end,:));
Ex_Alpha_ZZ = Ex_Alpha(:,9);
Ex_Alpha_Tr = sum(Ex_Alpha(:,[1,5,9]),2);

Sig_ZZZ     = Ex_Mu_Z.*Ex_Alpha_ZZ;
Sig_Norm    = Ex_Mu_Int.*Ex_Alpha_Tr;

% diaplay mode properties
Mode_List = [Ex_Ind,...
             Ex_Freq,...
             Sig_Norm,...
             Ex_Mu_Int,...
             Ex_Alpha_Tr,...
             Sig_ZZZ,...
             Ex_Mu_Z,...
             Ex_Alpha_ZZ,...
             ];

%% Update handles structure
handles.hMain          = GUI_Data_hModel.hMain;
handles.Structure      = Structure;
handles.MianGUI_Inputs = MainGUI_Inputs;
handles.OneDSFG        = OneDSFG;
guidata(hObject, handles);

% update the list on hPlot_Exciton GUI
set(handles.GUI_Modes.ModeList,'Data',Mode_List)

% call sorting to sort table with the same GUI setting
uitable_SortCallback(hObject, eventdata, handles)

function Update_Figure(hObject, eventdata, handles)
%% Re-assign variable names of GUI Inputs
GUI_handle  = handles.GUI_Modes;
Mode_Ind    = str2num(GUI_handle.Mode_Ind.String);
Mode_Type   = GUI_handle.Mode_Type.Value;
Plot_TDV    = GUI_handle.Plot_TDV.Value;
Scale_TDV   = str2double(GUI_handle.Scale_TDV.String);
Plot_Raman  = GUI_handle.Plot_Raman.Value;
Scale_Raman = str2double(GUI_handle.Scale_Raman.String);
Normalize   = GUI_handle.Normalize.Value;

%% Use Update Modes to update the structure and the corresponding Mu & Alpha 
Update_Modes(hObject, eventdata, handles)
handles = guidata(hObject);

Structure = handles.Structure;
OneDSFG   = handles.OneDSFG;

%% Calculate the exciton center
Center_Loc = Structure.center;

EigVecM   = OneDSFG.H.Sort_Ex_V;
EigVecM   = EigVecM(2:end,2:end).^2; % get ride of ground state
Center_Ex = EigVecM * Center_Loc;

%% Switch type of plotting mode

Num_Plot_Modes = length(Mode_Ind);

switch Mode_Type
    case 1
        Center = Center_Loc(Mode_Ind,:);
        Mu     = squeeze(OneDSFG.Mu.   Trans_Loc(1,Mode_Ind+1,:)); % shift by 1 to avoid ground state
        Alpha  = squeeze(OneDSFG.Alpha.Trans_Loc(1,Mode_Ind+1,:));
    case 2
        Center = Center_Ex(Mode_Ind,:);
        Mu     = squeeze(OneDSFG.Mu.   Trans_Ex(1,Mode_Ind+1,:)); % shift by 1 to avoid ground state
        Alpha  = squeeze(OneDSFG.Alpha.Trans_Ex(1,Mode_Ind+1,:));
end

% permute the matix dimension for spectial case
if Num_Plot_Modes == 1
    Mu    = Mu';
    Alpha = Alpha';
end

%% Plot molecule
if isgraphics(handles.hMain)
    GUI_Data_Main = guidata(handles.hMain);
    StructModel = get(GUI_Data_Main.GUI_Main.StructListBox,'Value');
else
    StructModel = 1;
end

[~,~,hPlotFunc] = StructureModel(StructModel);
hF = feval(hPlotFunc,Structure);
hAx = findobj(hF,'type','axes');
hold on

%% Plot Transition dipoles
if Plot_TDV
    if Normalize
        % normalize to unit vector for direction comparison
        Mu_Int = sqrt(sum(Mu.^2,2));
        Mu = bsxfun(@rdivide,Mu,Mu_Int);
    end
    Mu_S = Scale_TDV .* Mu; % Scale TDV vector in plot 
    
    quiver3(hAx,...
            Center(:,1),Center(:,2),Center(:,3),...
            Mu_S(:,1),Mu_S(:,2),Mu_S(:,3),0,...
            'LineWidth',2,...
            'Color',[255,128,0]./256);
end

%% plot Raman tensors
if Plot_Raman
    N_mesh   = 20;
    F_Color = [204,152,255]./255;
    if Normalize
        % normalize to unit vector for direction comparison
        Alpha_Tr = sum(Alpha(:,[1,5,9]),2);
        Alpha = bsxfun(@rdivide,Alpha,Alpha_Tr);
    end

    for i = 1: Num_Plot_Modes
        RamanM = reshape(Alpha(i,:),3,3);
        plot_Raman(hAx,RamanM,Center(i,:),Scale_Raman,N_mesh,F_Color)
    end
end
hold off

%% Figure setting
Fig_Title = ['Mode #: ' GUI_handle.Mode_Ind.String ];
hAx.Title.String = Fig_Title;

%% update handles
handles.Mode_Ind = Mode_Ind;
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
