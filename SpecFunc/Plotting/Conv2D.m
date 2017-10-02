function CVL = Conv2D(SG,GUI_Inputs)
% 
% This function convolute the input stick spectrum with spelected line
% shape.
% 

% Todo: integrate input parser with "Fig_Inputs"
%       Select which response to be plot
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.1  140723  Modified inputs
% 
% Ver. 1.0  140420  Isolated from "Plot2DSFG_DNA.m"
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

%% Debug
% FreqRange = 1550:1:1650;

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultFreqRange   = 1650:1750;
defaultLineShape   = 'Lorentzian';
defauleLineWidth   = 5;
defaultSpecType    = 'Absorptive';
defaultPathway     = 'All';

addOptional(INPUT,'FreqRange',defaultFreqRange);
addOptional(INPUT,'LineShape',defaultLineShape);
addOptional(INPUT,'LineWidth',defauleLineWidth);
addOptional(INPUT,'SpecType',defaultSpecType);
addOptional(INPUT,'Pathway',defaultPathway);

parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
FreqRange   = INPUT.Results.FreqRange;
LineShape   = INPUT.Results.LineShape;
LineWidth   = INPUT.Results.LineWidth;
Pathway     = INPUT.Results.Pathway;
SpecType    = INPUT.Results.SpecType;

%% Convert Sparse matrix back to full
MinF = FreqRange(1);
MaxF = FreqRange(end);

SG.SpecAccuR1  = full(SG.SpecAccuR1(MinF:MaxF,MinF:MaxF));
SG.SpecAccuR2  = full(SG.SpecAccuR2(MinF:MaxF,MinF:MaxF));
SG.SpecAccuR3  = full(SG.SpecAccuR3(MinF:MaxF,MinF:MaxF));
SG.SpecAccuNR1 = full(SG.SpecAccuNR1(MinF:MaxF,MinF:MaxF));
SG.SpecAccuNR2 = full(SG.SpecAccuNR2(MinF:MaxF,MinF:MaxF));
SG.SpecAccuNR3 = full(SG.SpecAccuNR3(MinF:MaxF,MinF:MaxF));

SG.Rephasing    = full(SG.Rephasing(MinF:MaxF,MinF:MaxF));
SG.NonRephasing = full(SG.NonRephasing(MinF:MaxF,MinF:MaxF));

%% FFT on selected pathway

switch Pathway
    case 'GB' % Ground state Bleach, R1
        Re_phasing_Res = SG.SpecAccuR1 ;
        NR_phasing_Res = SG.SpecAccuNR1;
        
    case 'SE' % Stimulated Emission, R2
        Re_phasing_Res = SG.SpecAccuR2 ;
        NR_phasing_Res = SG.SpecAccuNR2;
        
    case 'EA' % Excited state Absorption, R3
        Re_phasing_Res = -SG.SpecAccuR3 ;
        NR_phasing_Res = -SG.SpecAccuNR3;
        
    case 'All'
        Re_phasing_Res = SG.Rephasing   ;
        NR_phasing_Res = SG.NonRephasing;
        
end


%% Deal with Gaussian/Loentizan lineshape
% NumFreqPoint = numel(FreqRange);
% 
% center  = ceil(NumFreqPoint/2);
% [p1,p2] = meshgrid(1:NumFreqPoint,1:NumFreqPoint);
% 
% % Gaussian / Lorentizan line shape in frequency domain
% switch LineShape
%     case 'L'
%         FF = 1; % Counting for probe beam line width
%         lnshpf_R =((-1./(-(p2-center)+1i*LineWidth*FF)).*(1./((p1-center)+1i*LineWidth)));
%         lnshpf_N =((-1./( (p2-center)+1i*LineWidth*FF)).*(1./((p1-center)+1i*LineWidth)));
%     
%     case 'G'
%         lnshpf_R = ngaussval(sqrt((p1-center).^2+(p2-center).^2),LineWidth);
%         lnshpf_N = ngaussval(sqrt((p1-center).^2+(p2-center).^2),LineWidth);
%     case 'KK'
%         disp('not support KK lineshape in 2D yet...')
%     otherwise
%         lnshpf_R = 1;
%         lnshpf_N = 1;
%         disp('Plotting stick spectrum')     
% end
[ConvL,~] = Conv_LineShape(2,LineShape,FreqRange,LineWidth);
lnshpf_R  = ConvL.lnshpf_R;
lnshpf_N  = ConvL.lnshpf_N;

% export type of lineshape to output
CVL.Lineshape = LineShape;

%%
% Create Vis probe line shape
% -------------------------------------------------------------------------
% c = (3E8).*100./(1E12); % The speed of light (in cm/ps)
% delwn = 15; % This is the FWHM of the Lorentzian function (in cm^-1)
% sigma = c.*delwn./2; % This is the parameter sigma (in ps^-1)
% lortz = (sigma./pi).*(1./((sigma.^2) + ((p1-center).*c).^2)); % This is the Lorentzian function (in ps)
% normlortz = lortz./max(max(lortz)); % This is a normalized Lorentzian function (in arbitrary untis) to check that FWHM calculation is working correctly
% prb1 = normlortz;
% prb1t = ifftshift(ifft(fftshift(prb1),[],2));
% prbcor1 = (real(prb1t)./max(max(real(prb1t))));
% -------------------------------------------------------------------------        

% Convert the lineshape to time domain
% Rt  = ifftshift(ifft2(Re_phasing_Res));
% NRt = ifftshift(ifft2(NR_phasing_Res));
% lnshpt_R = ifftshift(ifft2(fftshift(lnshpf_R)));
% lnshpt_N = ifftshift(ifft2(fftshift(lnshpf_N)));
% CVL.R  = fft2(fftshift(lnshpt_R.*(Rt)));
% CVL.NR = fft2(fftshift(lnshpt_N.*(NRt)));

% frequency domain 2D convolution
CVL.R   = conv2(Re_phasing_Res,lnshpf_R,'same');
CVL.NR  = conv2(NR_phasing_Res,lnshpf_N,'same');
CVL.sum = CVL.R + CVL.NR;

% Not convoluted sticks
CVL.R_No_Conv   = Re_phasing_Res;
CVL.NR_No_Conv  = NR_phasing_Res;
CVL.sum_No_Conv = Re_phasing_Res + NR_phasing_Res;

% Output selected SpecType
switch SpecType
    case 'Absorptive'
        CVL.selected = CVL.sum;
        CVL.selected_No_Conv = CVL.sum_No_Conv;
    case 'Rephasing'
        CVL.selected = CVL.R;
        CVL.selected_No_Conv = CVL.R_No_Conv;
    case 'Non-rephasing'
        CVL.selected = CVL.NR;
        CVL.selected_No_Conv = CVL.NR_No_Conv;
end

% export the frequency array that the Specteal grid has
CVL.FreqRange = FreqRange;
