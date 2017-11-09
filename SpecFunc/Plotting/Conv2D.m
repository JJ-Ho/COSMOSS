function CVL = Conv2D(SSG,GUI_Inputs)

% This function convolute the input stick spectrum with spelected line
% shape.
% 

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
defaultFreqRange   = 1600:1750;
defaultPathway     = 'All';
defaultLineShape   = 'Lorentzian';
defauleLineWidth   = 5;
defaultSpecType    = 'Absorptive';

addOptional(INPUT,'FreqRange',defaultFreqRange);
addOptional(INPUT,'Pathway'  ,defaultPathway  );
addOptional(INPUT,'LineShape',defaultLineShape);
addOptional(INPUT,'LineWidth',defauleLineWidth);
addOptional(INPUT,'SpecType' ,defaultSpecType );
            
parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
FreqRange   = INPUT.Results.FreqRange;
Pathway     = INPUT.Results.Pathway;
LineShape   = INPUT.Results.LineShape;
LineWidth   = INPUT.Results.LineWidth;
SpecType    = INPUT.Results.SpecType;

%% Convert Sparse matrix back to full matrix
MinF = FreqRange(1);
MaxF = FreqRange(end);

SG.R1  = full(SSG.R1(MinF:MaxF,MinF:MaxF));
SG.R2  = full(SSG.R2(MinF:MaxF,MinF:MaxF));
SG.R3  = full(SSG.R3(MinF:MaxF,MinF:MaxF));
SG.NR1 = full(SSG.NR1(MinF:MaxF,MinF:MaxF));
SG.NR2 = full(SSG.NR2(MinF:MaxF,MinF:MaxF));
SG.NR3 = full(SSG.NR3(MinF:MaxF,MinF:MaxF));

SG.Rephasing    =  SG.R1 +  SG.R2 -  SG.R3;
SG.NonRephasing = SG.NR1 + SG.NR2 - SG.NR3;

%% Select pathways
switch Pathway
    case 'GB' % Ground state Bleach, R1
        Re_phasing_Res = SG.R1 ;
        NR_phasing_Res = SG.NR1;     
    case 'SE' % Stimulated Emission, R2
        Re_phasing_Res = SG.R2 ;
        NR_phasing_Res = SG.NR2;    
    case 'EA' % Excited state Absorption, R3
        Re_phasing_Res = -SG.R3 ;
        NR_phasing_Res = -SG.NR3;      
    case 'All'
        Re_phasing_Res = SG.Rephasing   ;
        NR_phasing_Res = SG.NonRephasing;    
end

%% Convolution
% Generate lineshapes
[ConvL2,~] = Conv_LineShape(2,LineShape,FreqRange,LineWidth);

% frequency domain 2D convolution
CVL.R   = conv2(Re_phasing_Res,ConvL2.lnshpf_R,'same');
CVL.NR  = conv2(NR_phasing_Res,ConvL2.lnshpf_N,'same');

% [ConvL1,~] = Conv_LineShape(1,'Lorentzian',FreqRange,20);
% ConvL1 = imag(ConvL1);
% N_Pump = size(Re_phasing_Res,1);
% for i = 1:N_Pump
%     CVL.R(i,:)  = conv(ConvL1,CVL.R(i,:),'same');
%     CVL.NR(i,:) = conv(ConvL1,CVL.NR(i,:),'same');
% end

CVL.sum = CVL.R + CVL.NR;

% Export Non-convoluted sticks
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

%% Other outputs
CVL.Lineshape = LineShape; % export type of lineshape to output
CVL.FreqRange = FreqRange; % export the frequency array
