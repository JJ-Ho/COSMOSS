function GUI = GUI_Plot_Modes(fig)
%% Pre-setting
Version = '1.4.2';

%% Add base layout
MainLayout = uix.VBoxFlex(...
    'Parent',fig,...
    'Spacing',1,...
    'Padding',3);

ModeListPanel = uix.BoxPanel( ...
    'Parent',MainLayout,...
    'Title','List of Modes and their properties','FontSize',14);

TDV_Raman_Panel = uix.BoxPanel( ...
    'Parent',MainLayout,...
    'Title','Transition dipole and Raman tensor','FontSize',14);

Response_Panel = uix.BoxPanel( ...
    'Parent',MainLayout,...
    'Title','Visulize molecular response','FontSize',14);

ButtonBox = uix.HBox('Parent',MainLayout,...
                     'Padding', 1, 'Spacing', 1);

set(MainLayout,'Heights',[-1,165,80,30])

%% ModeListPanel
SpecTypeList = {'FTIR','SFG','2DIR-2q','TwoDIR','TwoDSFG'};

ModeListLayout = uix.VBox('Parent',ModeListPanel,...
                           'Padding', 1, 'Spacing', 1);
% prep an empty table
uitable('Parent', ModeListLayout, ...
        'tag','ModeList',...
        'RearrangeableColumns','on',...
        'CellSelectionCallback',@(hObject,eventdata)Plot_Modes('uitable_CellSelection',hObject,eventdata,guidata(hObject)));

SpecTypeBox = uix.HBox('Parent',ModeListLayout,...
                       'Padding', 1, 'Spacing', 1);    
    uicontrol('Parent',SpecTypeBox,...
              'Style','text',...
              'String','Type: ',...
              'HorizontalAlignment','right',...
              'units', 'normalized',...
              'fontunits', 'point', 'fontsize', 14);
    uicontrol('Parent',SpecTypeBox,...
              'Style','popup',...
              'String',SpecTypeList,...
              'Value',2,...
              'tag','SpecType',...
              'HorizontalAlignment','left',...
              'units', 'normalized',...
              'fontunits', 'point', 'fontsize', 14);   
    uicontrol('Parent',SpecTypeBox,...
              'Style','popup',...
              'String',{},...
              'Value',1,...
              'tag','PathType',...
              'HorizontalAlignment','left',...
              'units', 'normalized',...
              'fontunits', 'point', 'fontsize', 14);   
    uicontrol('Parent', SpecTypeBox, ...
              'String', 'Refresh',...
              'fontunits', 'point', 'fontsize', 14,...
              'Callback',@(hObject,eventdata)Plot_Modes('Update_Modes',hObject,eventdata,guidata(hObject)));
    uix.Empty('Parent',SpecTypeBox);
    
    uicontrol('Parent',SpecTypeBox,...
              'Style','text',...
              'String','Sort by: ',...
              'HorizontalAlignment','right',...
              'units', 'normalized',...
              'fontunits', 'point', 'fontsize', 14);
    uicontrol('Parent',SpecTypeBox,...
              'Style','popup',...
              'String',{},...
              'Value',1,...
              'tag','SortInd',...
              'HorizontalAlignment','left',...
              'units', 'normalized',...
              'fontunits', 'point', 'fontsize', 14);
    uicontrol('Parent', SpecTypeBox, ...
              'String', 'Sort',...
              'fontunits', 'point', 'fontsize', 14,...
              'Callback',@(hObject,eventdata)Plot_Modes('uitable_Sort',hObject,eventdata,guidata(hObject)));

    uix.Empty('Parent',SpecTypeBox);
    
    set(SpecTypeBox,'Widths',[40,120,80,80,-1,55,120,60,-1])
    
set(ModeListLayout,'Heights',[-1,25])
         
%% TDV_Raman_Panel
PlotPanelLayout = uix.VBox('Parent',TDV_Raman_Panel,...
                           'Padding', 3, 'Spacing', 3);
                
    ModeBox = uix.HBox('Parent',PlotPanelLayout,...
                        'Padding', 1, 'Spacing', 1);
        uicontrol('Parent',ModeBox,...
                  'Style','checkbox',...
                  'String','Plot',...
                  'Value',1,...
                  'Tag','Mu_Alpha_Plot',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',ModeBox,...
                  'Style','popupmenu',...
                  'Tag','Mu_Alpha_Type',...
                  'Value',2,...
                  'String',{'Local mode', 'Exciton mode'},...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',ModeBox,...
                  'Style','text',...
                  'String',', index:',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',ModeBox,...
                  'Style','edit',...
                  'Tag','Mu_Alpha_Ind',...
                  'String','',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
    set(ModeBox,'Widths',[50,150,60,-1])
    
    TDV_ModeScaleBox = uix.HBox('Parent',PlotPanelLayout,...
                            'Padding', 1, 'Spacing', 1);
        uicontrol('Parent',TDV_ModeScaleBox,...
                  'Style','checkbox',...
                  'String','Transition Dipole Vectors,',...
                  'Value',1,...
                  'Tag','TDV_Plot',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14); 
        uicontrol('Parent',TDV_ModeScaleBox,...
                  'Style','text',...
                  'String','Scale:',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);     
        uicontrol('Parent',TDV_ModeScaleBox,...
                  'Style','edit',...
                  'String','1',...
                  'tag','TDV_Scale',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uix.Empty('Parent',TDV_ModeScaleBox);

    set(TDV_ModeScaleBox,'Widths',[200,60,40,-1])
        
    Raman_ModeScaleBox = uix.HBox('Parent',PlotPanelLayout,...
                            'Padding', 1, 'Spacing', 1);
        uicontrol('Parent',Raman_ModeScaleBox,...
                  'Style','checkbox',...
                  'String','Raman Tensors,',...
                  'Value',1,...
                  'Tag','Raman_Plot',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);     
        uicontrol('Parent',Raman_ModeScaleBox,...
                  'Style','text',...
                  'String','Scale:',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);     
        uicontrol('Parent',Raman_ModeScaleBox,...
                  'Style','edit',...
                  'String','1',...
                  'tag','Raman_Scale',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uix.Empty('Parent',Raman_ModeScaleBox);
        uicontrol('Parent',Raman_ModeScaleBox,...
                  'Style','text',...
                  'String','Display type:',...
                  'HorizontalAlignment','right',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);   
        uicontrol('Parent',Raman_ModeScaleBox,...
                  'Style','popupmenu',...
                  'Tag','Raman_Type',...
                  'Value',4,...
                  'String',{'Arrow3','Disk','Ellipsoid','Hyper Ellipsoid'},...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
    set(Raman_ModeScaleBox,'Widths',[200,60,40,-1,100,150])    

    uicontrol('Parent',PlotPanelLayout,...
              'Style','checkbox',...
              'String','Normalize modes for direction comparison,',...
              'Value',0,...
              'Tag','Normalize',...
              'HorizontalAlignment','left',...
              'units', 'normalized',...
              'fontunits', 'point', 'fontsize', 14);
                 
    EigenVecBox = uix.HBox('Parent',PlotPanelLayout,...
                        'Padding', 1, 'Spacing', 1);
        uicontrol('Parent',EigenVecBox,...
                  'Style','checkbox',...
                  'String','Plot eig. vec. of mode #:',...
                  'Value',0,...
                  'Tag','Plot_EigVec',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',EigenVecBox,...
                  'Style','edit',...
                  'Tag','EigVec_Ind',...
                  'String','',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',EigenVecBox,...
                  'Style','text',...
                  'String',', Scale:',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);     
        uicontrol('Parent',EigenVecBox,...
                  'Style','edit',...
                  'String','1',...
                  'tag','EigVec_Scale',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uix.Empty('Parent',EigenVecBox);
        uicontrol('Parent',EigenVecBox,...
                  'Style','text',...
                  'String','Convolute with:',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',EigenVecBox,...
                  'Style','popup',...
                  'String',{'None','Transition dipole','Raman Tensor'},...
                  'Value',1,...
                  'tag','EigVec_Conv',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);                                   

    set(EigenVecBox,'Widths',[180,50,60,40,-1,100,150])
          
set(PlotPanelLayout,'Heights',[-1,-1,-1,-1,-1])

%% Response_Panel
ResponseLayout = uix.VBox('Parent',Response_Panel,...
                           'Padding', 3, 'Spacing', 3);
    ThreeDBox = uix.HBox('Parent',ResponseLayout,...
                        'Padding', 1, 'Spacing', 1);
        uicontrol('Parent',ThreeDBox,...
                  'Style','checkbox',...
                  'String','Plot 3D map, ',...
                  'Value',1,...
                  'Tag','Sig_Plot3D',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);                    
        uicontrol('Parent',ThreeDBox,...
                  'Style','text',...
                  'String','Grid point: ',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',ThreeDBox,...
                  'Style','edit',...
                  'Tag','Sig_NGrid',...
                  'String','30',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',ThreeDBox,...
                  'Style','text',...
                  'String','Scale: ',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uicontrol('Parent',ThreeDBox,...
                  'Style','edit',...
                  'Tag','Sig_Scale',...
                  'String','1',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14);
        uix.Empty('Parent',ThreeDBox);
    set(ThreeDBox,'Widths',[150,80,40,60,40,-1])

    ContourBox = uix.HBox('Parent',ResponseLayout,...
                        'Padding', 1, 'Spacing', 1);
        uicontrol('Parent',ContourBox,...
                  'Style','checkbox',...
                  'String','Plot contour map',...
                  'Value',1,...
                  'Tag','Sig_PlotCT',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14); 
        uicontrol('Parent',ContourBox,...
                  'Style','checkbox',...
                  'String','Sum multiple modes',...
                  'Value',1,...
                  'Tag','Sig_PlotSum',...
                  'HorizontalAlignment','left',...
                  'units', 'normalized',...
                  'fontunits', 'point', 'fontsize', 14); 
        uix.Empty('Parent',ContourBox);
    set(ContourBox,'Widths',[150,180,-1])
    
    uix.Empty('Parent',ResponseLayout);
    
set(ResponseLayout,'Heights',[25,25,-1])

%% Button Box
uicontrol('Parent', ButtonBox, ...
          'String', 'Plot TDV/Raman',...
          'fontunits', 'point', 'fontsize', 14,...
          'Callback',@(hObject,eventdata)Plot_Modes('Update_TDV_Raman',hObject,eventdata,guidata(hObject)));

uicontrol('Parent', ButtonBox, ...
          'String', 'Plot Response',...
          'fontunits', 'point', 'fontsize', 14,...
          'Callback',@(hObject,eventdata)Plot_Modes('Update_Response',hObject,eventdata,guidata(hObject)));

uicontrol('Parent', ButtonBox, ...
          'String', 'Export handles',...
          'fontunits', 'point', 'fontsize', 14,...
          'Callback',@(hObject,eventdata)Plot_Modes('Export_Handle_Callback',hObject,eventdata,guidata(hObject)));
        
uicontrol('Parent',ButtonBox,...
          'Style','text',...
          'String',['Plot_Mode Version: v',Version],...
          'HorizontalAlignment','Center',...
          'units', 'normalized',...
          'fontunits', 'point', 'fontsize', 14);  
      
set(ButtonBox,'Width',[-1,-1,-1,200])
      
%% output handles
GUI = guihandles(fig);