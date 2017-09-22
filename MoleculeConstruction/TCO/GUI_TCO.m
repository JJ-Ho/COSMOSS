function GUI = GUI_TCO(fig) 
%% Pre-setting
Version = '1.2.0';

%% Add base layout
MainLayout = uix.VBoxFlex(...
    'Parent',fig,...
    'Spacing',1,...
    'Padding',5);

    OrientationPanel = uix.BoxPanel( ...
        'Parent',MainLayout,...
        'Title','Molecule Orientation','FontSize',14);

    LabelPanel = uix.BoxPanel( ...
        'Parent',MainLayout,...
        'Title','Local Mode Frequency','FontSize',14);

    FigurePanel = uix.BoxPanel( ...
        'Parent',MainLayout,...
        'Title','Molecule figure options','FontSize',14);
    
    ButtonBox = uix.VBoxFlex( ...
        'Parent', MainLayout,...
        'Padding', 2);

set(MainLayout,'Height',[150,120,50,-1])

%% Individule Chromophore Rotation 
OreintationLayout = uix.VBox('Parent',OrientationPanel,...
    'Padding', 5, 'Spacing', 5);

% oreintation inputs
Rot_HBox1 = uix.HBox(...
    'Parent',OreintationLayout,...
    'Spacing',3);

    uicontrol('Parent',Rot_HBox1,...
        'Style','text',...
        'String','Phi1');
    
    uicontrol('Parent',Rot_HBox1,...
        'Style','edit',...
        'String','0',...
        'Tag','Phi1',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

    uicontrol('Parent',Rot_HBox1,...
        'Style','text',...
        'String','Psi1');
    
    uicontrol('Parent',Rot_HBox1,...
        'Style','edit',...
        'String','0',...
        'Tag','Psi1',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

    uicontrol('Parent',Rot_HBox1,...
        'Style','text',...
        'String','Theta1');

    uicontrol('Parent',Rot_HBox1,...
        'Style','edit',...
        'String','0',...
        'Tag','Theta1',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));


Rot_HBox2 = uix.HBoxFlex(...
    'Parent',OreintationLayout,...
    'Spacing',3);

    uicontrol('Parent',Rot_HBox2,...
        'Style','text',...
        'String','Phi2');
    
    uicontrol('Parent',Rot_HBox2,...
        'Style','edit',...
        'String','0',...
        'Tag','Phi2',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

    uicontrol('Parent',Rot_HBox2,...
        'Style','text',...
        'String','Psi2');
    
    uicontrol('Parent',Rot_HBox2,...
        'Style','edit',...
        'String','0',...
        'Tag','Psi2',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

    uicontrol('Parent',Rot_HBox2,...
        'Style','text',...
        'String','Theta2');

    uicontrol('Parent',Rot_HBox2,...
        'Style','edit',...
        'String','0',...
        'Tag','Theta2',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

%% Translation inputs for the second chromophore
TransVecLayout = uix.HBox('Parent',OreintationLayout,...
    'Padding', 5, 'Spacing', 5);

    uicontrol('Parent',TransVecLayout,...
        'Style','text',...
        'String','Trans Vec.');

    uicontrol('Parent',TransVecLayout,...
        'Style','edit',...
        'String','5,0,0',...
        'Tag','Trans',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

set(TransVecLayout,'Widths',[80,-1])

%% Rotation of the whole two-chromophore system

% oreintation inputs
RotAll_HBox1 = uix.HBox(...
    'Parent',OreintationLayout,...
    'Spacing',3);

    uicontrol('Parent',RotAll_HBox1,...
        'Style','text',...
        'String','Rot. All');

    uicontrol('Parent',RotAll_HBox1,...
        'Style','text',...
        'String','X');
    
    uicontrol('Parent',RotAll_HBox1,...
        'Style','edit',...
        'String','0',...
        'Tag','Rot_X',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

    uicontrol('Parent',RotAll_HBox1,...
        'Style','text',...
        'String','Y');
    
    uicontrol('Parent',RotAll_HBox1,...
        'Style','edit',...
        'String','0',...
        'Tag','Rot_Y',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

    uicontrol('Parent',RotAll_HBox1,...
        'Style','text',...
        'String','Z');

    uicontrol('Parent',RotAll_HBox1,...
        'Style','edit',...
        'String','0',...
        'Tag','Rot_Z',...
        'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));
    
set(RotAll_HBox1,'Widths',[20,-1,-1,-1,-1,-1,-1])

%% Labeling panel
LabelPanelLayout = uix.HBoxFlex('Parent',LabelPanel,...
    'Padding', 5, 'Spacing', 5);

    Freq_TextBox = uix.VBox('Parent',LabelPanelLayout);
        uicontrol('Parent',Freq_TextBox,...
        'Style','text',...
        'String','Non-Labeled Freq:',...
        'HorizontalAlignment','right',...
        'units', 'normalized',...
        'fontunits', 'point', 'fontsize', 14);
    
        uicontrol('Parent',Freq_TextBox,...
        'Style','text',...
        'String','Labeled Freq:',...
        'HorizontalAlignment','right',...
        'units', 'normalized',...
        'fontunits', 'point', 'fontsize', 14);
    
        uicontrol('Parent',Freq_TextBox,...
        'Style','text',...
        'String','Anharmonicty:',...
        'HorizontalAlignment','right',...
        'units', 'normalized',...
        'fontunits', 'point', 'fontsize', 14);
    
        uicontrol('Parent',Freq_TextBox,...
        'Style','text',...
        'String','Index of modes:',...
        'HorizontalAlignment','right',...
        'units', 'normalized',...
        'fontunits', 'point', 'fontsize', 14);
    
    

    Freq_EditBox = uix.VBox('Parent',LabelPanelLayout);
        uicontrol('Parent',Freq_EditBox,...
        'Style','edit',...
        'String','1650',...
        'Tag','NLFreq',...
        'units', 'normalized',...
        'fontunits', 'point', 'fontsize', 14);
    
        uicontrol('Parent',Freq_EditBox,...
        'Style','edit',...
        'String','1640',...
        'Tag','LFreq',...
        'units', 'normalized',...
        'fontunits', 'point', 'fontsize', 14);
    
        uicontrol('Parent',Freq_EditBox,...
        'Style','edit',...
        'String','10',...
        'Tag','Anharm',...
        'units', 'normalized',...
        'fontunits', 'point', 'fontsize', 14);   
        
        uicontrol('Parent',Freq_EditBox,...
        'Style','edit',...
        'String','None',...
        'Tag','L_Index',...
        'units', 'normalized',...
        'fontunits', 'point', 'fontsize', 14);

set(LabelPanelLayout,'Widths',[130,-1])    
   
%% Plot molecule options
FigurePanelLayout = uix.HBox('Parent',FigurePanel,...
    'Padding', 5, 'Spacing', 5);

        uicontrol('Parent',FigurePanelLayout,...
                  'Style','checkbox',...
                  'String','Atoms,',...
                  'Value',1,...
                  'Tag','Plot_Atoms',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);     
        uicontrol('Parent',FigurePanelLayout,...
                  'Style','checkbox',...
                  'String','Bonds.',...
                  'Value',1,...
                  'Tag','Plot_Bonds',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',FigurePanelLayout,...
                  'Style','checkbox',...
                  'String','Axis.',...
                  'Value',1,...
                  'Tag','Plot_Axis',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uix.Empty('Parent',FigurePanelLayout);
    set(FigurePanelLayout,'Widths',[80,80,80,-1])

%% Button
uicontrol( 'Parent', ButtonBox, ...
'String', 'Generate',...
'Callback',@(hObject,eventdata)Model_TCO('UpdateStructure',hObject,eventdata,guidata(hObject)));

uicontrol( 'Parent', ButtonBox, ...
'String', 'Plot Molecule',...
'Callback',@(hObject,eventdata)Model_TCO('PlotMolecule',hObject,eventdata,guidata(hObject)));

uicontrol( 'Parent', ButtonBox, ...
'String', 'Plot Modes',...
'Callback',@(hObject,eventdata)Model_TCO('PlotModes',hObject,eventdata,guidata(hObject)));

uicontrol( 'Parent', ButtonBox, ...
'String', 'Export handles',...
'Callback',@(hObject,eventdata)Model_TCO('Export_Handle_Callback',hObject,eventdata,guidata(hObject)));

uicontrol('Parent',ButtonBox,...
          'Style','text',...
          'String',['Model Version: v',Version],...
          'HorizontalAlignment','Center',...
          'units', 'normalized',...
          'fontunits', 'point', 'fontsize', 14);       
      
set(ButtonBox,'Height',[-1,-1,-1,-1,20])

%% output handles
GUI = guihandles(fig);

    