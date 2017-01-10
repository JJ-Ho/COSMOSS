function hF = H_Cut(Cut_F,GUI_Inputs,TwoD_Data)
%% Debug Prep varaibles
% Cut_F = 1707;
% GUI_Inputs = ParseGUI_Main(Data_COSMOSS.hGUIs);
% 
% % TwoD_Data = Data_COSMOSS.TwoDSFG;
% TwoD_Data = Data_COSMOSS.TwoDIR;

%% Extract info from 2D data
FreqRange  = GUI_Inputs.FreqRange;
PlotCursor = GUI_Inputs.PlotCursor;

C  = TwoD_Data.CVL.selected;
SS = TwoD_Data.CVL.selected_No_Conv;
SI = TwoD_Data.Int;
F  = TwoD_Data.Freq;

%% Extract convoluted H-Cut
C_Cut_X = FreqRange;
F_Min = FreqRange(1);

C_Ind = Cut_F - F_Min + 1;
C_Cut_Y = real(C(C_Ind,:));

%% Extract individule Sticks
Type = fieldnames(SI);

SI_M_P = zeros(length(C_Cut_X),7); % first column zero for accumulation
SI_M_N = zeros(length(C_Cut_X),7);

for i = 1:length(Type)
    Tmp_F = abs(F.(Type{i}));
    Tmp_SI_Ind = eq(Tmp_F(:,1),Cut_F);
    
    Tmp_SI_X = Tmp_F(Tmp_SI_Ind,3);
    
    if or(eq(i,3),eq(i,6))
        Tmp_Y = -SI.(Type{i})(Tmp_SI_Ind); % for R3/NR3
    else
        Tmp_Y = SI.(Type{i})(Tmp_SI_Ind);
    end
    
    P_Ind = ge(Tmp_Y,0);
    N_Ind = lt(Tmp_Y,0);
    
    X_P_Sub = Tmp_SI_X(P_Ind)-F_Min+1;
    X_N_Sub = Tmp_SI_X(N_Ind)-F_Min+1;
    
    Y_P_Sub = Tmp_Y(P_Ind);
    Y_N_Sub = Tmp_Y(N_Ind);
    
    % remove terms out of range
    NG = length(C_Cut_X);
    P_remove_Ind = or(X_P_Sub>NG,X_P_Sub<=0);
    N_remove_Ind = or(X_N_Sub>NG,X_N_Sub<=0);
    
    X_P_Sub(P_remove_Ind) = [];
    Y_P_Sub(P_remove_Ind) = [];
    X_N_Sub(N_remove_Ind) = [];
    Y_N_Sub(N_remove_Ind) = [];
    
    if isempty(X_P_Sub)
        SI_M_P(:,i+1) = zeros(length(C_Cut_X),1);
    else
        SI_M_P(:,i+1) = accumarray(X_P_Sub,Y_P_Sub,[NG,1]);
    end
    
    if isempty(X_N_Sub)
        SI_M_N(:,i+1) = zeros(length(C_Cut_X),1);
    else
        SI_M_N(:,i+1) = accumarray(X_N_Sub,Y_N_Sub,[NG,1]);
    end
    
end

%% Extract summed sticks
SS_X = FreqRange;
SS_Y = SS(C_Ind,:);

SS_Ind = gt(abs(SS_Y),0);

SS_X = SS_X(SS_Ind);
SS_Y = SS_Y(SS_Ind);

%% Draw figure
% plot convoluted curve
hF = figure;

hAx_SI = subplot(3,1,[1,2]);
hold(hAx_SI,'on')
    
    % plot convolution curve
    yyaxis(hAx_SI,'left')
    plot(hAx_SI,C_Cut_X,C_Cut_Y)

    % plot baseline for individule sticks
    line(hAx_SI,[C_Cut_X(1);C_Cut_X(end)],[0;0],'Color',[0,0,0]);

    % Plot individule sticks
    LineWidth = 4;
    LColor = [255,  0,  0;...     % R1
                0,255,  0;...     % R2
                0,  0,255;...     % R3
              255,  0,255;...     % NR1
                0,153,  0;...     % NR2
                0,255,255]./255;  % NR3

    yyaxis(hAx_SI,'right')
    for j = 1:length(Type)    
        Tmp_SI_X   = [C_Cut_X',C_Cut_X'];
        Tmp_SI_Y_P = [sum(SI_M_P(:,1:j),2),sum(SI_M_P(:,1:j+1),2)];
        Tmp_SI_Y_N = [sum(SI_M_N(:,1:j),2),sum(SI_M_N(:,1:j+1),2)];

        % Set intensity cut of to speed up 
        Int_CutOff = 0.05;
        CutOff_I_P = gt(abs(SI_M_P(:,j+1)),max(abs(SI_M_P(:,j+1)))*Int_CutOff);
        CutOff_I_N = gt(abs(SI_M_N(:,j+1)),max(abs(SI_M_N(:,j+1)))*Int_CutOff);

        line(hAx_SI,...
                Tmp_SI_X(CutOff_I_P,:)',...
             Tmp_SI_Y_P(CutOff_I_P,:)',...
            'Color',LColor(j,:),'LineWidth',LineWidth,'LineStyle','-');

        line(hAx_SI,...
               Tmp_SI_X(CutOff_I_N,:)',...
             Tmp_SI_Y_N(CutOff_I_N,:)',...
            'Color',LColor(j,:),'LineWidth',LineWidth,'LineStyle','-');

    end
hold(hAx_SI,'off')

hAx_SS = subplot(3,1,3);
hold(hAx_SS,'on')

    yyaxis(hAx_SS,'left')
    % plot convolution curve
    plot(hAx_SS,C_Cut_X,C_Cut_Y)
    
    % plot baseline for individule sticks
    line(hAx_SS,[C_Cut_X(1);C_Cut_X(end)],[0;0],'Color',[1,0,0]);

    % Draw stick sum
    Tmp_SS_X = [SS_X',SS_X'];
    Tmp_SS_Y = [zeros(length(SS_Y),1),SS_Y'];
    % Set intensity cut of to speed up 
    CutOff_I_SS = gt(abs(SS_Y(:)),max(abs(SS_Y(:)))*Int_CutOff);
    
    yyaxis(hAx_SS,'right')
    line(hAx_SS,...
         Tmp_SS_X(CutOff_I_SS,:)',...
         Tmp_SS_Y(CutOff_I_SS,:)',...
        'Color',[0,0,0],'LineWidth',LineWidth,'LineStyle','-');

hold(hAx_SS,'off')

%% Figure adjustment
hF.Units = 'normalized'; % use normalized scale

% Individule Stick subplot
hAx_SI.FontSize = 14;
hAx_SI.XLim = [C_Cut_X(1);C_Cut_X(end)];
hAx_SI.XTickLabel = '';
hAx_SI.XGrid = 'on';
hAx_SI.YGrid = 'on';
hAx_SI.XMinorGrid = 'on';

yyaxis(hAx_SI,'left')
YLim_max = 1.1*max(abs(hAx_SI.YLim));
hAx_SI.YLim = [-YLim_max;YLim_max];
hAx_SI.YLabel.String = 'Convoltion Intensity';

yyaxis(hAx_SI,'right')
YLim_max = 1.1*max(abs(hAx_SI.YLim));
hAx_SI.YLim = [-YLim_max;YLim_max];
hAx_SI.YLabel.String = 'Stick Intensity';

% Stick sum subplot
hAx_SS.FontSize = 14;
hAx_SS.XLim = [C_Cut_X(1);C_Cut_X(end)];
hAx_SS.XLabel.String = 'cm^{-1}';
hAx_SS.XGrid = 'on';
hAx_SS.YGrid = 'on';
hAx_SS.XMinorGrid = 'on';

yyaxis(hAx_SS,'left')
YLim_max = 1.1*max(abs(hAx_SS.YLim));
hAx_SS.YLim = [-YLim_max;YLim_max];

yyaxis(hAx_SS,'right')
YLim_max = 1.1*max(abs(hAx_SS.YLim));
hAx_SS.YLim = [-YLim_max;YLim_max];

% deal with title
SpecType = TwoD_Data.SpecType;
Title_String = [SpecType,' H-Cut @ ',num2str(Cut_F),' cm^{-1}'];

if PlotCursor    
    % Call pointer
    SI.fh = hF;
    SI.ax = hAx_SI;
    Pointer_N(SI) % use normalized scale
    
    title(hAx_SS,Title_String,'FontSize',16);
else
    title(hAx_SI,Title_String,'FontSize',16);
end

