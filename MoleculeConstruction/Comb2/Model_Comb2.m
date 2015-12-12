function varargout = Model_Comb2(varargin)
% MODEL_COMB2 MATLAB code for Model_Comb2.fig
%      MODEL_COMB2, by itself, creates a new MODEL_COMB2 or raises the existing
%      singleton*.
%
%      H = MODEL_COMB2 returns the handle to a new MODEL_COMB2 or the handle to
%      the existing singleton*.
%
%      MODEL_COMB2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODEL_COMB2.M with the given input arguments.
%
%      MODEL_COMB2('Property','Value',...) creates a new MODEL_COMB2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Model_Comb2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Model_Comb2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Model_Comb2

% Last Modified by GUIDE v2.5 12-Dec-2015 12:22:30

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


% --- Executes just before Model_Comb2 is made visible.
function Model_Comb2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Model_Comb2 (see VARARGIN)

% Choose default command line output for Model_Comb2
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
    disp('Running in stand alone mode.')    
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Model_Comb2 wait for user response (see UIRESUME)
% uiwait(handles.Model_Comb2);


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
[hStructure1,~,hPlotFunc1] = StructureModel(StructModel,handles);

handles.hStruc1    = hStructure1;
handles.hPlotFunc1 = hPlotFunc1;
guidata(hObject,handles)

function Struc2(hObject, eventdata, handles)

GUI_Struc = handles.GUI_Struc;
StructModel    = get(GUI_Struc.StructBox2,'Value');
[hStructure2,~,hPlotFunc2] = StructureModel(StructModel,handles);

handles.hStruc2    = hStructure2;
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
Rot_Phi   = GUI_Inputs.Rot_Phi;
Rot_Psi   = GUI_Inputs.Rot_Psi;
Rot_Theta = GUI_Inputs.Rot_Theta;

TransV = [Trans_X,Trans_Y,Trans_Z];
RM = R1_ZYZ_0(Rot_Phi,Rot_Psi,Rot_Theta);

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
Comb2.center       = M_Center;
Comb2.freq         = M_Freq;
Comb2.anharm       = M_Anharm;
Comb2.mu           = M_TDV;
Comb2.alpha        = M_Raman;
Comb2.alpha_matrix = M_Raman_Matrix;
% Comb2.AtomSerNo    = M_AtomSerNo;
Comb2.Num_Modes    = size(M_TDV,1);
Comb2.XYZ          = M_XYZ;
Comb2.FilesName    = 'Comb2';

Structure = Comb2;
%% export back to handles and Main GUI if any
Data_Main.Structure = Structure;

handles.StrucData1  = StrucData1;
handles.StrucData2  = StrucData2;
handles.Comb2       = Comb2;
handles.Data_Main   = Data_Main;

if isfield(handles,'hMain')
    guidata(handles.hMain,Data_Main)
end
guidata(hObject,handles)

function PlotComb2_Callback(hObject, eventdata, handles)

%% retreive data from handles
hStruc1        = handles.hStruc1;
hStruc2        = handles.hStruc2;
StrucData1     = handles.StrucData1;
StrucData2     = handles.StrucData2;

hPlotFunc1     = handles.hPlotFunc1;
hPlotFunc2     = handles.hPlotFunc2;

%% Retreive GUI inputs
GUI_Struc = handles.GUI_Struc;
GUI_Inputs = ParseGUI_Comb2(GUI_Struc);
Trans_X   = GUI_Inputs.Trans_X;
Trans_Y   = GUI_Inputs.Trans_Y;
Trans_Z   = GUI_Inputs.Trans_Z;
Rot_Phi   = GUI_Inputs.Rot_Phi;
Rot_Psi   = GUI_Inputs.Rot_Psi;
Rot_Theta = GUI_Inputs.Rot_Theta;

TransV = [Trans_X,Trans_Y,Trans_Z];
RM = R1_ZYZ_0(Rot_Phi,Rot_Psi,Rot_Theta);

Center2 = StrucData2.center;
COM2 = sum(Center2,1)./size(Center2,1);

%% make figure
hF1 = feval(hPlotFunc1,StrucData1);
hF2 = feval(hPlotFunc2,StrucData2);
hAx1 = findobj(hF1,'type','axes');
hAx2 = findobj(hF2,'type','axes');

hFcomb2 = figure;
hAx_Comb2 = axes;
copyobj(allchild(hAx1),hAx_Comb2);
copyobj(allchild(hAx2),hAx_Comb2);

close(hF1)
close(hF2)
hold on 
% plot the XYZ axis of the second structure
Scale = 2;
Axis2_0 = [1,0,0;0,1,0;0,0,1];
Axis2_R = RM*Axis2_0 .* Scale;

Origin   = bsxfun(@times,COM2,[1,1,1]');

CMatrix = [1,0,0;0,1,0;0,0,1];

for i = 1:3
quiver3(Origin(i,1),Origin(i,2),Origin(i,3),...
        Axis2_R(i,1),Axis2_R(i,2),Axis2_R(i,3),0,...
        'LineWidth',5,...
        'Color',CMatrix(i,:));
end  
    

hold off
%% figure options
axis equal;
rotate3d on

xlabel('X')
ylabel('Y')

function Export_Handle_Callback(hObject, eventdata, handles)
% export handles back to work space
assignin('base', 'H', handles)
disp('Updated handles exported!')
