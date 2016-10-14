function Output = PlotFTIR(PDB_Data,GUI_Inputs)
%% PlotFTIR
%  
% Given one exciton alpha, mu matrix and respective one exciton state
% energy (cm-1), this function can generate 1DSFG plots with different
% polarization combinations {xx,yy,zz,xy,yz,xz}*{x,y,z} = 18 plots.
% 
% 
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.4  140922  Add Output part;
%                   Add Inputparser
% 
% Ver. 1.3  140717  Add Frequency axis GUI read in part
% 
% Ver. 1.2  140616  Add CouplingOpt
% 
% Ver. 1.1  140605  Isolate from PlotOneDSFG.m
% 
% Ver. 1.0  130729  Isolated from TwoDSFG_Simulation.
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

%% Debug
% PDB_Data = GetAcid;
% handles.PDB_Data = PDB_Data;
% 
% GUI_Inputs.debug = 1;
% GUI_Inputs.Coupling = 'NN';

%% Inputs parser
% Turn Output from Read GUI to cell array
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

% Default values
% defaultLabel_Index  = 'Non';
% defaultLabel_Freq   = 1700;
% defaultPlotStick    = 1;
defaultCouplingType = 'TDC';
defaultBeta_NN      = 0.8;
% defaultF_Min        = 1600;
% defaultF_Max        = 1800;
% defaultLineWidth    = 5;
% defaultPlotCursor   = 0;
% defaultLineShape    = 'G';

% add Optional inputs / Parameters
% addOptional(INPUT,'Label_Index' ,defaultLabel_Index);
% addOptional(INPUT,'Label_Freq'  ,defaultLabel_Freq);
% addOptional(INPUT,'PlotStick'   ,defaultPlotStick);
addOptional(INPUT,'CouplingType',defaultCouplingType);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);
% addOptional(INPUT,'F_Min'       ,defaultF_Min);
% addOptional(INPUT,'F_Max'       ,defaultF_Max);
% addOptional(INPUT,'LineWidth'   ,defaultLineWidth);
% addOptional(INPUT,'PlotCursor'  ,defaultPlotCursor);
% addOptional(INPUT,'LineShape'   ,defaultLineShape);

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
% PlotStick    = INPUT.Results.PlotStick;
CouplingType = INPUT.Results.CouplingType;
Beta_NN      = INPUT.Results.Beta_NN;
% F_Min        = INPUT.Results.F_Min;
% F_Max        = INPUT.Results.F_Max;
% LineWidth    = INPUT.Results.LineWidth;
% PlotCursor   = INPUT.Results.PlotCursor;
% LineShape    = INPUT.Results.LineShape;

%% Main
Num_Modes = PDB_Data.Num_Modes;

H = ExcitonH(PDB_Data,'ExMode','OneEx','CouplingType',CouplingType,'Beta_NN',Beta_NN);

mu = MuAlphaGen(PDB_Data,H,'Mode','Mu');

muEx = mu.Trans_Ex;
muEx_Vec = muEx(2:Num_Modes+1,1,:);

Freq01Ex = H.Sort_Ex_Freq;
freq_OneD = Freq01Ex(2:Num_Modes+1);

mu2_OneD = sum(muEx_Vec.^2,3); % E-field of FTIR signal is mu^2 base on feynmann duagram! 

%% Generate figure 
% mu_OneD = mu2_OneD;
% PlotStick    = 1;
% F_Min        = 1550;
% F_Max        = 1700;
% LineWidth    = 5;
% PlotCursor   = 0;
% LineShape    = 'L';
% 
% 
% hF = figure; hold on
% 
% if eq(PlotStick,1)
%     line([freq_OneD';freq_OneD'],[zeros(1,Num_Modes);mu_OneD'])
% end
% 
% % Get Frequency axis range
% spec_range = F_Min:F_Max;
% 
% spec_array1 = bsxfun(@times,ones(Num_Modes,length(spec_range)),spec_range);
% spec_array2 = bsxfun(@minus,spec_array1,freq_OneD);
% 
% 
% switch LineShape 
%     case 'G' % Gaussian
%         LineShape = exp(-(spec_array2.^2)./(LineWidth^2));
%         CVL = bsxfun(@times,LineShape,mu_OneD); 
%         CVL_Total = sum(CVL,1);
%     case 'L' % Lorentzain 
%         LineWidth = LineWidth/2;
%         LineShape = LineWidth./((spec_array2.^2)+(LineWidth^2));
%         CVL = bsxfun(@times,LineShape,mu_OneD); 
%         CVL_Total = sum(CVL,1); 
%     case 'KK'
%         disp('KK is not support for FTIR...')
%         CVL_Total = zeros(size(spec_range));
% end
% 
% 
% % Normalize the convoluted lineshape to maximum stick height.
% Norm = max(abs(mu_OneD(:)))./max(abs(CVL_Total));
% CVL_Total = CVL_Total.*Norm;
% 
% plot(spec_range,CVL_Total,'-')
% hold off
% 
% %% figure setting 
% hF.Units = 'normalized'; % use normalized scale
% hAx = hF.CurrentAxes;
% hAx.FontSize = 14;
% hAx.XLim = [spec_range(1),spec_range(end)];
% hAx.XLabel.String = 'cm^{-1}';
% 
% grid on
% 
% if PlotCursor
%     % Call pointer
%     S.fh = hF;
%     S.ax = hAx;
%     Pointer_N(S) % use normalized scale
% else
%     FilesName     = PDB_Data.FilesName;
%     FilesName_Reg = regexprep(FilesName,'\_','\\_');
%     Coupling_Reg  = regexprep(CouplingType,'\_','\\_');
%     Title_String  = ['FTIR, ',FilesName_Reg,', Coupling:',Coupling_Reg];
%     title(Title_String,'FontSize',16);    
% end
% 
% % % integrate the curve area
% % Area = trapz(spec_range,Gaussian_Toatl);
% % StickSum = sum(mu_OneD);
% % 
% % Title = [num2str(Area), ',  ', num2str(StickSum)];
% % title(Title,'FontSize',16);

%% Output
Output.Num_Modes    = Num_Modes;
Output.Response1D   = mu2_OneD;
Output.freq_OneD    = freq_OneD;
Output.SpecType     = 'FTIR';
Output.FilesName    = PDB_Data.FilesName;
Output.CouplingType = CouplingType;


Output.H  = H;
Output.Mu = mu;
% Output.hF = hF;
% Output.X  = spec_range;
% Output.Y  = CVL_Total;
