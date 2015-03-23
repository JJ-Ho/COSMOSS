function Plot2DIR(CVL,FreqRange,H,Mu_Ex)
% 
% This function plot 2DIR and other information

% ------- Version log -----------------------------------------------------
%  
% Ver. 1.1  141014  Use Pointer_N instead of Pointer for normailized unit
% 
% Ver. 1.0  140723  Isolated from "TwoDIR_Main.m"
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014

f = figure('Position',[600,200,500,450],...
           'units','pix');

% Generate Mu value table and eigen mode frequencies
Mod_Ine_All = 2:3; % for two coupled oscillators

F = H.Sort_Ex_Freq(Mod_Ine_All)';


Mu_Vec = squeeze(Mu_Ex(1,Mod_Ine_All,:));
Mu_Vec_L = sqrt(sum(Mu_Vec.^2,2))';
T = [ Mu_Vec_L.^2;
      Mu_Vec_L.^4;
      Mu_Vec(:,1)';
      Mu_Vec(:,2)';
      Mu_Vec(:,3)'];

Beta = ones(1,2).*H.Beta(1,2);
  
dat = [F;Beta;T]'; 
rnames = {'S1','S2'};
cnames = {'Freq','Beta','IR_Int','2DIR_Int','Mu_X','Mu_Y','Mu_Z'};

t = uitable('Parent',f,...
            'Data',dat,...
            'ColumnFormat',{'bank','bank','bank','bank','bank','bank','bank'},...
            'ColumnWidth' ,{60,50,60,80,40,40,40},...
            'ColumnName',cnames,... 
            'RowName',rnames,...
            'units','pix',...
            'fontsize',13,...
            'Position',[40,10,420,60]);
        
% Make 2DSFG plot
h1 = axes('units','pix','Position',[100,80,350,350]);
% h1 = axes('units','pix','Position',[100,100,300,300]);
set(gcf,'CurrentAxes',h1)

% Num_Contour = 20;
% [C,Ch] = contour(FreqRange,FreqRange,real(CVL.selected),Num_Contour,'LineWidth',2);
% shading flat
% set(gca,'DataAspectRatio',[1 1 1])

Z = -1:(2/40):1;
% V = Z;
V = [Z(1:20) Z(22:41)];

% contour(FreqRange,FreqRange,-1.*real(CVL.selected)./max(max(abs(real(CVL.selected)))),V,'LineWidth',1.5)
contour(FreqRange,FreqRange,-1.*real(CVL.selected),40,'LineWidth',1.5)
% caxis([-1 1]);
shading flat
set(gca,'DataAspectRatio',[1 1 1])

% Plot diagonal line
hold on; plot(FreqRange,FreqRange); hold off

% Set colorbar
colorbar
Amp = max(abs(caxis));
caxis([-Amp Amp])

% Call pointer
S.fh = f;
S.ax = h1;
% Pointer(S)
Pointer_N(S)

% figure; plot(FreqRange,diag(real(CVL.selected)))