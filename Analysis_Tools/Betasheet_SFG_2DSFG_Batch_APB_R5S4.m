% batch script for 1DSFG/2DSFG generation of ideal betasheet
%% Inputs
%- Major options  ---------------------------------------------------------
Scan_Theta = 0:10:350;
% Scan_Theta = 10;
O.Coupling = 'Cho_APB';

%- Figure options ---------------------------------------------------------
Fig_Save = 0;
saveGIF  = 1;
BaseFileName = 'APB_R5S4_';
PathName = pwd;

Struc_XLim = [-12.00, 12.00];
Struc_YLim = [-11.40, 11.40];
Struc_ZLim = [-15.00, 15.00];
Struc_View = [-120,5];

OneDSFG_YLim      = [-0.6,0.6];
TwoDSFG_CLim      = [-1E3,1E3];
TwoDSFG_Diag_CLim = [-6E4,6E4];

Struc_XLimMode        = 'Manual';
Struc_YLimMode        = 'Manual';
Struc_ZLimMode        = 'Manual';
OneDSFG_YLimMode      = 'Manual';
TwoDSFG_CLimMode      = 'Manual';
TwoDSFG_Diag_CLimMode = 'Manual';

Mode_Ind         = [4,10,6,2];
Mode_Stick_Width = 5;

Scale_TDV   = 20;
Scale_Raman = 0.5;
N_mesh      = 20;

%- Structural part --------------------------------------------------------
StrucInput.SheetType = 2; % 'Anti'
StrucInput.N_Residue = 5;
StrucInput.N_Strand  = 4;
StrucInput.Trans_X   = 0;
StrucInput.Trans_Y   = 0;
StrucInput.Trans_Z   = 4.75;
StrucInput.Twist_X   = 0;
StrucInput.Twist_Y   = 0;
StrucInput.Twist_Z   = 0;
StrucInput.Phi_D     = 0;
StrucInput.Psi_D     = 0;
StrucInput.Theta_D   = 0;
StrucInput.NLFreq    = 1644;
StrucInput.LFreq     = 1604;
StrucInput.Anharm    = 20;
StrucInput.L_Index   = 'None';

% For 1D ------------------------------------------------------------------
O.P_IR         = 0;
O.P_Vis1D      = 0;
O.P_Sig1D      = 0;
O.A_IR         = 90;
O.A_Vis1D      = 90;
O.A_Sig1D      = 90;
% ------------------------------------------------------------------------

% For 2D ------------------------------------------------------------------
O.A_Pump       = 90;
O.A_Probe      = 90;
O.A_Vis2D      = 90;
O.A_Sig2D      = 90;
O.P_Pump1      = 0;
O.P_Pump2      = 0;
O.P_Probe      = 0;
O.P_Vis2D      = 0;
O.P_Sig2D      = 0;
% ------------------------------------------------------------------------
  
O.Avg_Rot      = 1; %'Phi' C_Inf
O.Avg_Mirror   = 1; % no mirror plane

O.Num_Contour  = 20;
O.PlotStick    = 1;
O.PlotNorm     = 1;
O.PlotCursor   = 0;

O.F_Min     = 1550;
O.F_Max     = 1700;
O.FreqRange = O.F_Min:O.F_Max;

O.LineShape = 'L';
O.SpecType = 'Abs';
O.Pathway = 'All';
O.Signal_Type = 'Hetero';


%% Scanning 
 

for i = 1:length(Scan_Theta)

    StrucInput.Theta_D = Scan_Theta(i);
    %% Structure 
    BB        = ConstuctBetaSheet(StrucInput);
    Structure = GetAmideI(BB.Num_Atoms,...
                          BB.XYZ,...
                          BB.AtomName,...
                          BB.FilesName,...
                          StrucInput);

    Structure.N_Residue         = StrucInput.N_Residue;
    Structure.N_Strand          = StrucInput.N_Strand;
    Structure.N_Mode_per_Starnd = StrucInput.N_Residue-1;

    % C terminus Index
    Structure.Ind_H = BB.Ind_H;
    Structure.Ind_O = BB.Ind_O;

    % Betasheet orientation info export
    Structure.TransV = [StrucInput.Trans_X,StrucInput.Trans_Y,StrucInput.Trans_Z];
    Structure.TwistV = [StrucInput.Twist_X,StrucInput.Twist_Y,StrucInput.Twist_Z];
    Structure.RotV   = [StrucInput.Phi_D  ,StrucInput.Psi_D  ,StrucInput.Theta_D];

    hF_Struc = Plot_Betasheet_AmideI(Structure);
    hAx_Struc = findobj(hF_Struc,'type','axes');
    
    hAx_Struc.XLim = Struc_XLim;
    hAx_Struc.YLim = Struc_YLim;
    hAx_Struc.ZLim = Struc_ZLim;
    hAx_Struc.XLimMode = Struc_XLimMode;
    hAx_Struc.YLimMode = Struc_YLimMode;
    hAx_Struc.ZLimMode = Struc_ZLimMode;
    
    view(Struc_View)
    
    %% OneDSFG
    OneDSFG    = OneDSFG_Main(Structure,O);
    O.PlotNorm = 0;
    hF_OneDSFG = PlotOneDSFG(OneDSFG,O);
    hAx_OneDSFG = findobj(hF_OneDSFG,'type','axes');
    
    % Apply Color and Broader linewith to selected modes
    hSticks = findobj(hAx_OneDSFG,'Tag','Stick');
    Num_Plot_Modes = length(Mode_Ind);
    UnWanted = [0,0,0;1,1,1;1,0,0;0,1,0;0,0,1];
    Mode_colors = distinguishable_colors(Num_Plot_Modes,UnWanted);
    
    for m = 1:Num_Plot_Modes
        Stick_Ind = Structure.Num_Modes +1 - Mode_Ind(m);
        hSticks(Stick_Ind).Color = Mode_colors(m,:);
        hSticks(Stick_Ind).LineWidth = Mode_Stick_Width;
    end
    
    hAx_OneDSFG.YLim = OneDSFG_YLim;
    hAx_OneDSFG.YLimMode = OneDSFG_YLimMode;
    hAx_OneDSFG.XGrid = 'on';
    %hAx_OneDSFG.YGrid = 'on';
    
    %% FTIR
    FTIR = PlotFTIR(Structure,O);
    FTIR_X = FTIR.X;
    FTIR_Y = FTIR.Y;
    FTIR_Y_Norm = FTIR_Y./abs(max(FTIR_Y));
    Scale_OneDSFG = max(abs(hAx_OneDSFG.YLim));
    FTIR_Y_Scale = FTIR_Y_Norm.*Scale_OneDSFG;
    
    hold(hAx_OneDSFG,'on')
    plot(hAx_OneDSFG,FTIR_X,FTIR_Y_Scale,'r','LineWidth',2)
    hold(hAx_OneDSFG,'off')
    
    %% Plot selected Exiton Modes
    hold(hAx_Struc,'on')
    
    % Calculate the exciton center
    Center_Loc = Structure.center;

    EigVecM   = OneDSFG.H.Sort_Ex_V;
    EigVecM   = EigVecM(2:end,2:end).^2; % get ride of ground state
    Center_Ex = EigVecM * Center_Loc;
    
    Center = Center_Ex(Mode_Ind,:);
    Mu     = squeeze(OneDSFG.Mu.   Trans_Ex(1,Mode_Ind+1,:)); % shift by 1 to avoid ground state
    Alpha  = squeeze(OneDSFG.Alpha.Trans_Ex(1,Mode_Ind+1,:));
    
    Mu_S = Scale_TDV .* Mu; % Scale TDV vector in plot 
    
    for r = 1: Num_Plot_Modes
        
        quiver3(hAx_Struc,...
            Center(r,1),Center(r,2),Center(r,3),...
            Mu_S(r,1),Mu_S(r,2),Mu_S(r,3),0,...
            'LineWidth',Mode_Stick_Width,...
            'Color',Mode_colors(r,:));
        
        RamanM = reshape(Alpha(r,:),3,3);
        plot_Raman(hAx_Struc,RamanM,Center(r,:),Scale_Raman,N_mesh,Mode_colors(r,:))
        
    end
        
    hold(hAx_Struc,'off')
    
    %% Two DSFG
    [SpectraGrid,Response] = TwoDSFG_Main(Structure,O);

    CVL = Conv2D(SpectraGrid,O);
    CVL.FilesName = Structure.FilesName; % pass filesname for figure title
    hF_TwoDSFG = Plot2DSFG(CVL,O);
    hAx_TwoDSFG = findobj(hF_TwoDSFG,'type','axes');
    hAx_TwoDSFG_cbar = findobj(hF_TwoDSFG,'type','colorbar');
    
    if strcmp(TwoDSFG_CLimMode,'Manual')
        hAx_TwoDSFG.CLim = TwoDSFG_CLim;
        hAx_TwoDSFG.CLimMode = TwoDSFG_CLimMode;
    end
    
    hAx_TwoDSFG.XGrid ='on';
    hAx_TwoDSFG.YGrid ='on';
    
    %% TwoDSFG diagonal sticks
    Re_Diag = diag(SpectraGrid.Rephasing);
    NR_Diag = diag(SpectraGrid.NonRephasing);
    
    hF_2D_Diag = figure;
    line([O.FreqRange(1);O.FreqRange(end)],[0;0],'Color',[0,0,0])
    line([O.FreqRange;O.FreqRange],[zeros(1,length(O.FreqRange));(Re_Diag+NR_Diag)'],'LineWidth',4)
    
    hAx_2D_Diag = findobj(hF_2D_Diag,'type','axes');
    
    if strcmp(TwoDSFG_Diag_CLimMode,'Manual')
        line([O.FreqRange(end);O.FreqRange(end)],TwoDSFG_Diag_CLim,'Color',[1,0,0],'LineWidth',4)
        hAx_2D_Diag.YLim = TwoDSFG_Diag_CLim;
        hAx_2D_Diag.YLimMode = TwoDSFG_Diag_CLimMode;
    else 
        line([O.FreqRange(end);O.FreqRange(end)],hAx_2D_Diag.YLim,'Color',[1,0,0],'LineWidth',4)
    end
    
    hAx_2D_Diag.Color = 'none';
    hAx_2D_Diag.Visible = 'off';
    
    % Plot Rephasing / Nonrephasing separate
    %line([O.FreqRange;O.FreqRange],[zeros(1,length(O.FreqRange));Re_Diag'],'LineWidth',2)
    %line([O.FreqRange;O.FreqRange],[zeros(1,length(O.FreqRange));NR_Diag'],'LineWidth',2)
    
    % plot in log scale
    %Log_Re_Diag = log(abs(Re_Diag)).*sign(Re_Diag);
    %Log_NR_Diag = log(abs(Re_Diag)).*sign(NR_Diag);
    %Line([O.FreqRange;O.FreqRange],[zeros(1,length(O.FreqRange));Log_Re_Diag'],'LineWidth',2)
    %line([O.FreqRange;O.FreqRange],[zeros(1,length(O.FreqRange));Log_NR_Diag'],'LineWidth',2)
    

    %% Paste 2DSFG diag sticks to 2DSFG figure
    hAx_2D_Diag_sub = copyobj(hAx_2D_Diag,hF_TwoDSFG);
    
    pos_2D = hAx_TwoDSFG.Position;
    pos_2D_Diag = pos_2D;
    pos_2D_Diag(1) = pos_2D_Diag(1)+0.037;
    pos_2D_Diag(3) = pos_2D_Diag(4)*0.750;
    pos_2D_Diag(4) = pos_2D_Diag(4)*0.300;
    set(hAx_2D_Diag_sub,'Position',pos_2D_Diag);
    
    %% Makeing subplot
    hF = figure;
    
    hAx_Struc_sub   = copyobj(hAx_Struc,hF);
    hAx_OneDSFG_sub = copyobj(hAx_OneDSFG,hF);
    hAx_TwoDSFG_sub = copyobj([hAx_TwoDSFG,hAx_TwoDSFG_cbar],hF);
    hAx_2D_Diag_sub = copyobj(hAx_2D_Diag,hF);
    
    pos_Struct  = [0.0343    0.1068    0.2614    0.7864];
    pos_1DSFG   = [0.3320    0.1259    0.2760    0.7626];
    pos_2DSFG   = [0.6733    0.1025    0.2503    0.8150];
    pos_cBar    = [0.9466    0.1083    0.0149    0.7929];
    pos_Base    = [0.0690    0.4886    0.7726    0.4010];
    pos_2D_Diag = [0.6788    0.1070    0.2393    0.3150];
    
    set(hAx_Struc_sub     ,'Position',pos_Struct);
    set(hAx_OneDSFG_sub   ,'Position',pos_1DSFG);
    set(hAx_TwoDSFG_sub(1),'Position',pos_2DSFG);
    set(hAx_TwoDSFG_sub(2),'Position',pos_cBar);
    set(hAx_2D_Diag_sub   ,'Position',pos_2D_Diag);
    
    colormap('jet')

    set(hF,'Unit','Normalized','Position',pos_Base)
    
    Frame_all(i) = getframe(hF);
    
    %% Save figure
    SavePath = PathName;
    if Fig_Save
        
        SaveNamePrefix = [BaseFileName,'Th',num2str(Scan_Theta(i))];
        FigName1  = [SaveNamePrefix,'_Structure'];
        SaveName1 = [SavePath,'/', FigName1];
        FigName2  = [SaveNamePrefix,'_SFG_FTIR'];
        SaveName2 = [SavePath,'/', FigName2];
        FigName3  = [SaveNamePrefix,'_2DSFG'];
        SaveName3 = [SavePath,'/', FigName3];
        FigName4  = [SaveNamePrefix,'_2DSFG_Diag'];
        SaveName4 = [SavePath,'/', FigName4];
        
        saveas(hF_Struc  ,SaveName1,'png') 
        saveas(hF_OneDSFG,SaveName2,'png')
        saveas(hF_TwoDSFG,SaveName3,'png')
        %saveas(hF_2D_Diag,SaveName4,'png')

    end
    
    if saveGIF
        close all
    end
end

disp('Loop finished...')

if saveGIF
    %% Save Frame
    SaveName = [SavePath,'/',BaseFileName];
    FrameName = [SaveName,'Frame_0-360'];
    save(FrameName,'Frame_all') 

    %% save gif
    load(FrameName)
    filename = [SaveName,'0-360.gif'];
    DT = 0.5;

    for j = 1:length(Scan_Theta)

        frame = Frame_all(j);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);

        if j == 1;
          imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',DT);
        else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DT);
        end

    end
    disp('GIF saved...')
end
    %% Plot frame
% load('Frame_all_0-360.mat')
% hF = figure; set(hF,'Position',[207 494 2051 527]);set(gca,'Position',[0,0,1,1]),movie(Frame_all,5,1)

% load(FrameName)
% hF = figure; set(hF,'Unit','Normalized','Position',[0.0690    0.4886    0.7726    0.4010]);set(gca,'Position',[0,0,1,1]),movie(Frame_all,5,1)
