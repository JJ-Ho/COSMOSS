function varargout = Model_Comb2(varargin)
% HMODEL MATLAB code for hModel.fig
%      HMODEL, by itself, creates a new HMODEL or raises the existing
%      singleton*.
%
%      H = HMODEL returns the handle to a new HMODEL or the handle to
%      the existing singleton*.
%
%      HMODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HMODEL.M with the given input arguments.
%
%      HMODEL('Property','Value',...) creates a new HMODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_Comb2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_Comb2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hModel

% Last Modified by GUIDE v2.5 28-Feb-2016 23:01:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Model_Comb2_OpeningFcn, ...
                   'gui_OutputFcn',  @Model_Comb2_OutputFcn, ...
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

% --- Executes just before hModel is made visible.
function Model_Comb2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hModel (see VARARGIN)

% Choose default command line output for hModel
handles.output = hObject;

% Call createInterface to create GUI elements
GUI_Struc = GUI_Comb2(hObject);
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
    disp('Running Model_Comb2 in stand alone mode.')    
end
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Model_Comb2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Struc1(hObject, eventdata, handles)

GUI_Struc = handles.GUI_Struc;
StructModel    = get(GUI_Struc.StructBox1,'Value');
[hStructure1,~,hPlotFunc1] = StructureModel(StructModel);

hStruc1 = feval(hStructure1,handles.hMain);

handles.hStruc1    = hStruc1;
handles.hPlotFunc1 = hPlotFunc1;
guidata(hObject,handles)

function Struc2(hObject, eventdata, handles)

GUI_Struc = handles.GUI_Struc;
StructModel    = get(GUI_Struc.StructBox2,'Value');
[hStructure2,~,hPlotFunc2] = StructureModel(StructModel);

hStruc2 = feval(hStructure2,handles.hMain);

handles.hStruc2    = hStruc2;
handles.hPlotFunc2 = hPlotFunc2;
guidata(hObject,handles)

function UpdateStructure(hObject, eventdata, handles)
GUI_Struc = handles.GUI_Struc;

hStruc1        = handles.hStruc1;
hStruc2        = handles.hStruc2;
StrucGUI_Data1 = guidata(hStruc1);
StrucGUI_Data2 = guidata(hStruc2);
StrucData1     = StrucGUI_Data1.Structure;
StrucData2     = StrucGUI_Data2.Structure;

Data_Main = handles.Data_Main;
%% Retreive GUI inputs
GUI_Inputs = ParseGUI_Comb2(GUI_Struc);
Trans_X   = GUI_Inputs.Trans_X;
Trans_Y   = GUI_Inputs.Trans_Y;
Trans_Z   = GUI_Inputs.Trans_Z;
Rot_Phi   = GUI_Inputs.Rot_Phi/180*pi;
Rot_Psi   = GUI_Inputs.Rot_Psi/180*pi;
Rot_Theta = GUI_Inputs.Rot_Theta/180*pi;

TransV = [Trans_X,Trans_Y,Trans_Z];
% RM = R1_ZYZ_0(Rot_Phi,Rot_Psi,Rot_Theta);
RM = Rx(Rot_Phi)*Ry(Rot_Psi)*Rz(Rot_Theta);

%% Shift the center of mass of each structure to origin
Center1 = StrucData1.center;
COM1 = sum(Center1,1)./size(Center1,1);

Center2 = StrucData2.center;
COM2 = sum(Center2,1)./size(Center2,1);

%% Move the second molecule and merge
% Center 
Center1_0 = bsxfun(@minus,Center1,COM1);
Center2_0 = bsxfun(@minus,Center2,COM2);

XYZ2_R = (RM*Center2_0')';
Center2_T = bsxfun(@plus,XYZ2_R,TransV);

M_Center = [Center1_0; Center2_T];

% mu 
TDV1 = StrucData1.mu;
TDV2 = StrucData2.mu;

TDV2_R = (RM*TDV2')';

M_TDV = [TDV1;TDV2_R];

% alpha matrix
Raman_Matrix1 = StrucData1.alpha_matrix;
Raman_Matrix2 = StrucData2.alpha_matrix;
Num_Modes2     = StrucData2.Num_Modes;

Raman_Matrix2_R = zeros(size(Raman_Matrix2));
for i=1:Num_Modes2
    Raman_Matrix2_R(i,:,:) = RM*squeeze(Raman_Matrix2(i,:,:))*RM';
end

M_Raman_Matrix = cat(1,Raman_Matrix1,Raman_Matrix2_R);

% alpha vector
Raman1   = reshape(Raman_Matrix1  ,[],9);
Raman2_R = reshape(Raman_Matrix2_R,[],9);
M_Raman  = reshape(M_Raman_Matrix ,[],9);

% XYZ 
XYZ1 = StrucData1.XYZ;
XYZ2 = StrucData2.XYZ;

XYZ1_0 = bsxfun(@minus,XYZ1,COM1);
XYZ2_0 = bsxfun(@minus,XYZ2,COM2);

XYZ2_R = (RM*XYZ2_0')';
XYZ2_T = bsxfun(@plus,XYZ2_R,TransV);

M_XYZ = [XYZ1_0;XYZ2_T];

%% update Structure1 and 2

StrucData1.center       = Center1_0;
StrucData1.mu           = TDV1;
StrucData1.alpha_matrix = Raman_Matrix1;
StrucData1.alpha        = Raman1;
StrucData1.XYZ          = XYZ1_0;

StrucData2.center       = Center2_T;
StrucData2.mu           = TDV2_R;
StrucData2.alpha_matrix = Raman_Matrix2_R;
StrucData2.alpha        = Raman2_R;
StrucData2.XYZ          = XYZ2_T;

%% Merge each catagory that doesn't need to move
% freq
Freq1 = StrucData1.freq;
Freq2 = StrucData2.freq;
M_Freq = [Freq1;Freq2];

% anharmonicity
Anharm1 = StrucData1.anharm;
Anharm2 = StrucData2.anharm;
M_Anharm = [Anharm1;Anharm2];

% % Atom Serial Number
% EAtomSerNo  = StrucData1.AtomSerNo;
% Shift_Index = size(EXYZ,1);
% FAtomSerNo  = StrucData2.AtomSerNo + Shift_Index;
% M_AtomSerNo = [EAtomSerNo;FAtomSerNo];

%% Output
Structure.center       = M_Center;
Structure.freq         = M_Freq;
Structure.anharm       = M_Anharm;
Structure.mu           = M_TDV;
Structure.alpha        = M_Raman;
Structure.alpha_matrix = M_Raman_Matrix;
Structure.Num_Modes    = size(M_TDV,1);
Structure.XYZ          = M_XYZ;
Structure.FilesName    = 'Comb2';

% Export into Structure so it can be passsed around different GUIs
Structure.StrucData1  = StrucData1;
Structure.StrucData2  = StrucData2;
Structure.StructModel = 5;

%% export back to handles and Main GUI if any
handles.Structure   = Structure;

Data_Main.Structure = Structure;
handles.Data_Main   = Data_Main;

if isfield(handles,'hMain')
    guidata(handles.hMain,Data_Main)
        
    % change Name of Main GUI to help identifying which Structural Model is
    % using
    Model_Name    = handles.hModel.Name;
    handles.hMain.Name = ['COSMOSS: ' Model_Name];
end
guidata(hObject,handles)

function hF = PlotMolecule(hObject, eventdata, handles)
GUI_Struc = handles.GUI_Struc;
GUI_Inputs = ParseGUI_Comb2(GUI_Struc);

hF = PlotComb2(handles,GUI_Inputs);

UpdateStructure(hObject, eventdata, handles)

function PlotModes(hObject, eventdata, handles)
Plot_Modes(handles.hModel);

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'hModel_Comb2', handles)
disp('Updated handles exported!')
