%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
            set(hGUIs_COSMOSS.F_Min  ,'String','1550')
            set(hGUIs_COSMOSS.F_Max  ,'String','1750')
            
            GUI_data.hCOSMOSS = hCOSMOSS;
            
            disp('Running Model_PDB_AmideI directly from COSMOSS...')
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

% Prep necessary data to be saved in GUI_data
GUI_data.hModel_PDB_AmideI = hModel_PDB_AmideI;
GUI_data.hGUIs             = hGUIs;

% Update handles structure
guidata(hModel_PDB_AmideI, GUI_data);

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
%^ GUI Setup ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



function LoadStructure(hObject, eventdata, GUI_data)
%% Read GUI variables
hGUIs  = GUI_data.hGUIs;
COSMOSS_GUI_Inputs = ParseGUI_AmideI(hGUIs);
Preprocessed = COSMOSS_GUI_Inputs.Preprocessed;

%% Get pdb file location 
PWD = pwd;
PDB_Path = [PWD, '/StructureFiles/PDB/'];

[FilesName,PathName,~] = uigetfile({'*.pdb','PDB file'; ...
                                    '*,*','All Files'},...
                                    'MultiSelect','on',...
                                    'Select inputs',PDB_Path);
                                
%% Parse molecule structure
if Preprocessed
    % read the pre-processed MD sanpshots with the same molecule
    TextPattern = 'ATOM %*f %s %*s %*s %*f %f %f %f %*f %*f %*s';
    if iscell(FilesName)
        N_File = size(FilesName,2);
        fid = fopen([PathName FilesName{1}]);
        ATest = textscan(fid,TextPattern,'CollectOutput',1);
        fclose(fid);
        XYZ = zeros([size(ATest{2}),N_File]); % Creat all zeros RawData Matrix
        
        progressbar;
        for i=1:N_File
            progressbar(i/N_File)
            % read preprocessed pdb file
            fid = fopen([PathName FilesName{i}]);
            A = textscan(fid,TextPattern,'CollectOutput',1);
            fclose(fid);

            XYZ(:,:,i) = A{2};
        end
        AtomName  = A{1};
        Num_Atoms = size(XYZ,1);
        
        C = strsplit(PathName,'/');
        FilesName  = C{end-1};
        
        % deal with the GUI inputs in COSMOSS
        hCOSMOSS = GUI_data.hCOSMOSS;
        Data_COSMOSS = guidata(hCOSMOSS);
        COSMOSS_hGUIs = Data_COSMOSS.hGUIs;
        COSMOSS_hGUIs.Sampling.Value = 1;
        COSMOSS_hGUIs.Sample_Num.String = num2str(N_File);
        COSMOSS_hGUIs.FWHM.String = num2str(0);
        
    else
        N_File = 1;
        fid = fopen([PathName FilesName]);
        A = textscan(fid,TextPattern,'CollectOutput',1);
        fclose(fid);
        
        AtomName   = A{1};
        XYZ(:,:,1) = A{2};
        Num_Atoms = size(XYZ,1);
    end
else 
    N_File = 1;
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
end

%% output to GUI
GUI_data.PDB.Num_Atoms = Num_Atoms;
GUI_data.PDB.XYZ       = XYZ;
GUI_data.PDB.AtomName  = AtomName;
GUI_data.PDB.FilesName = FilesName;
GUI_data.PDB.N_File    = N_File;

guidata(hObject,GUI_data)

%% Update PDB name and data to GUI
set(GUI_data.hGUIs.PDB_Name,'String',FilesName)
UpdateStructure(hObject, eventdata, GUI_data)

function UpdateStructure(hObject, eventdata, GUI_data)
%% Read GUI variables
hGUIs  = GUI_data.hGUIs;
GUI_Inputs = ParseGUI_AmideI(hGUIs);
    
%% Construct molecule
Num_Atoms = GUI_data.PDB.Num_Atoms;
XYZ       = GUI_data.PDB.XYZ;
AtomName  = GUI_data.PDB.AtomName;
FilesName = GUI_data.PDB.FilesName;
N_File    = GUI_data.PDB.N_File;

% test # modes and pre-allocate matix
Tmp1 = GetAmideI(XYZ(:,:,1),AtomName,FilesName,GUI_Inputs);
Nmodes = Tmp1.Nmodes;
Tmp_LocMu     = zeros(Nmodes,3,N_File);
Tmp_LocAlpha  = zeros(Nmodes,9,N_File);
Tmp_LocCenter = zeros(Nmodes,3,N_File);
Tmp_XYZ       = zeros(Num_Atoms,3,N_File);

for i = 1:N_File
    Tmp = GetAmideI(XYZ(:,:,i),AtomName,FilesName,GUI_Inputs);
    Tmp_LocMu(:,:,i)     = Tmp.LocMu;
    Tmp_LocAlpha(:,:,i)  = Tmp.LocAlpha;
    Tmp_LocCenter(:,:,i) = Tmp.LocCenter;
    Tmp_XYZ(:,:,i)       = Tmp.XYZ;
end

Structure = StructureData;
Structure.XYZ       = Tmp_XYZ;
Structure.AtomName  = Tmp1.AtomName;
Structure.COM       = Tmp1.COM;

Structure.LocCenter = Tmp_LocCenter;
Structure.LocFreq   = Tmp1.LocFreq;
Structure.LocAnharm = Tmp1.LocAnharm;
Structure.LocMu     = Tmp_LocMu;
Structure.LocAlpha  = Tmp_LocAlpha;

Structure.FilesName = Tmp1.FilesName;
Structure.Extra.AmideIAtomSerNo = Tmp1.Extra.AmideIAtomSerNo;

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

hF = PlotXYZfiles_AmideI(GUI_data.Structure,GUI_Inputs);

function PlotModes(hObject, eventdata, GUI_data)
Plot_Modes(GUI_data.hModel_PDB_AmideI);






