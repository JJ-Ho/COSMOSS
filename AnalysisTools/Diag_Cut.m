function hF = Diag_Cut(GUI_Inputs,TwoD_Data)
%% Debug Prep varaibles
% Cut_F = 1707;
% 
% GI = ParseGUI_Main(Data_COSMOSS.hGUIs);
% FreqRange = GI.FreqRange;
% 
% % TwoD_Data = Data_COSMOSS.TwoDSFG;
% TwoD_Data = Data_COSMOSS.TwoDIR;

%% Extract info from 2D data
FreqRange  = GUI_Inputs.FreqRange;
PlotCursor = GUI_Inputs.PlotCursor;

C = TwoD_Data.CVL.selected;
S = TwoD_Data.Int;
F = TwoD_Data.Freq;

%% Extract convoluted H-Cut
C_Cut_X = FreqRange;
F_Min = FreqRange(1);

C_Cut_Y = real(diag(C));

%% Extract individule Sticks
Type = fieldnames(S);

S_M_P = zeros(length(C_Cut_X),7); % first column zero for accumulation
S_M_N = zeros(length(C_Cut_X),7);

for i = 1:length(Type)
    Tmp_F = abs(F.(Type{i}));
    Tmp_S_Ind = eq(Tmp_F(:,1),Cut_F);
    
    Tmp_X = Tmp_F(Tmp_S_Ind,3);
    
    if or(eq(i,3),eq(i,6))
        Tmp_Y = -S.(Type{i})(Tmp_S_Ind); % for R3/NR3
    else
        Tmp_Y = S.(Type{i})(Tmp_S_Ind);
    end
    
    P_Ind = ge(Tmp_Y,0);
    N_Ind = lt(Tmp_Y,0);
    
    X_P_Sub = Tmp_X(P_Ind)-F_Min+1;
    X_N_Sub = Tmp_X(N_Ind)-F_Min+1;
    
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
        S_M_P(:,i+1) = zeros(length(C_Cut_X),1);
    else
        S_M_P(:,i+1) = accumarray(X_P_Sub,Y_P_Sub,[NG,1]);
    end
    
    if isempty(X_N_Sub)
        S_M_N(:,i+1) = zeros(length(C_Cut_X),1);
    else
        S_M_N(:,i+1) = accumarray(X_N_Sub,Y_N_Sub,[NG,1]);
    end
    
end

%% Extract summed sticks

%% Draw figure
% plot convoluted curve
hF = figure;
yyaxis left
plot(C_Cut_X,C_Cut_Y)
hold on 

% plot baseline
line([C_Cut_X(1);C_Cut_X(end)],[0;0],'Color',[0,0,0]);

% Plot sticks
LineWidth = 4;
LColor = [255,  0,  0;...     % R1
            0,255,  0;...     % R2
            0,  0,255;...     % R3
          255,  0,255;...     % NR1
            0,153,  0;...     % NR2
            0,255,255]./255;  % NR3

% LScale = max(C_Cut_Y)/max([sum(S_M_P,2);sum(-S_M_N,2)]);
LScale = 1;
yyaxis right
for j = 1:length(Type)    
    Tmp_X   = [C_Cut_X',C_Cut_X'];
    Tmp_Y_P = [sum(S_M_P(:,1:j),2),sum(S_M_P(:,1:j+1),2)].*LScale;
    Tmp_Y_N = [sum(S_M_N(:,1:j),2),sum(S_M_N(:,1:j+1),2)].*LScale;
   
    % Set intensity cut of to speed up 
    CutOff_V = 0.01;
    CutOff_I_P = gt(abs(S_M_P(:,j+1)),max(abs(S_M_P(:,j+1)))*CutOff_V);
    CutOff_I_N = gt(abs(S_M_N(:,j+1)),max(abs(S_M_N(:,j+1)))*CutOff_V);

    line(  Tmp_X(CutOff_I_P,:)',...
         Tmp_Y_P(CutOff_I_P,:)',...
        'Color',LColor(j,:),'LineWidth',LineWidth,'LineStyle','-');
    
    line(  Tmp_X(CutOff_I_N,:)',...
         Tmp_Y_N(CutOff_I_N,:)',...
        'Color',LColor(j,:),'LineWidth',LineWidth,'LineStyle','-');

end
hold off

%% Figure adjustment
hAx = findobj(hF,'Type','Axes');

hAx.FontSize = 14;
hAx.XLim = [C_Cut_X(1);C_Cut_X(end)];
hAx.XLabel.String = 'cm^{-1}';
hAx.XGrid = 'on';
hAx.YGrid = 'on';
hAx.XMinorGrid = 'on';

yyaxis left
YLim_max = 1.1*max(abs(hAx.YLim));
hAx.YLim = [-YLim_max;YLim_max];
hAx.YLabel.String = 'Convoltion Intensity';

yyaxis right
YLim_max = 1.1*max(abs(hAx.YLim));
hAx.YLim = [-YLim_max;YLim_max];
hAx.YLabel.String = 'Stick Intensity';

% deal with title
hF.Units = 'normalized'; % use normalized scale
SpecType = TwoD_Data.SpecType;
Title_String{1} = [SpecType,' H-Cut @ ',num2str(Cut_F),' cm^{-1}'];

if PlotCursor
    Title_String{2} = '';
    % Call pointer
    S.fh = hF;
    S.ax = hAx;
    Pointer_N(S) % use normalized scale
end

title(Title_String,'FontSize',16);

