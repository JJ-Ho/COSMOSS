function Output = Bin2D_Vec(Response,FreqRange,SignalColumn)
%% Bin2DSFG
% 
% Given the response, slected polarization column (response value) and the
% binning range (FreqRange), this script bin the input data into a square
% matrix with 1cm-1 bin size and the matrix is of 
% (lenght of FreqRange) * (lenght of FreqRange) size.
% 
% ------- Version log -----------------------------------------------------
% 
% Ver. 2.1  130911  Redefine R1-NR3 to Response.Bin(R1-NR3) to unifiy 
%                   different calculation need. Such as TwoDSFG_Thymine and
%                   TwoDSFG_AmideI.
%                   Rename from plot2DSFG ro Bin2DSFG
% 
% Ver. 2.0  130910  Redefine R1-NR3 as Lab frame response with <R> and J
% 
% Ver. 1.3  130901  Move cutoff to "TwoDSFG_Simulation" to accommadate
%                   Rotational avgerage. Also reassign response to Averaged
%                   Response functions.
% 
% Ver. 1.2  130806  Move cutoff to "Feynman_2DSFG_ForLoopV2" to reduce
%                   output text size.
% 
% Ver. 1.1  130805  Functionalized 
% 
% Ver. 1.0  130805  This is a vectorized version that come from fruitful
%                   discussion with Yu, Kuang and Shi, Liang.
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013


R1 = Response.BinR1;
R2 = Response.BinR2;
R3 = Response.BinR3;
NR1 = Response.BinNR1;
NR2 = Response.BinNR2;
NR3 = Response.BinNR3;

% ReverseFreq = FreqRange(end:-1:1);
% NumFreq = numel(FreqRange);

%% size and polarization

NumResR1 = size(R1,1);
NumResR2 = size(R2,1);
NumResR3 = size(R3,1);
NumResNR1 = size(NR1,1);
NumResNR2 = size(NR2,1);
NumResNR3 = size(NR3,1);

ResponseR1 = R1(:,SignalColumn);
ResponseR2 = R2(:,SignalColumn);
ResponseR3 = R3(:,SignalColumn);
ResponseNR1 = NR1(:,SignalColumn);
ResponseNR2 = NR2(:,SignalColumn);
ResponseNR3 = NR3(:,SignalColumn);

%% round to 1cm-1 resolution
pumpR1 = round(R1(:,1));
probR1 = round(R1(:,3));
pumpR2 = round(R2(:,1));
probR2 = round(R2(:,3));
pumpR3 = round(R3(:,1));
probR3 = round(R3(:,3));

pumpNR1 = round(NR1(:,1));
probNR1 = round(NR1(:,3));
pumpNR2 = round(NR2(:,1));
probNR2 = round(NR2(:,3));
pumpNR3 = round(NR3(:,1));
probNR3 = round(NR3(:,3));

%% R1
[~,PumpR1Ref]= ndgrid(ones(NumResR1,1),FreqRange);
[~,ProbR1Ref]= ndgrid(ones(NumResR1,1),FreqRange);

PumpR1MinusRef = bsxfun(@plus,PumpR1Ref,pumpR1)'; % since the pump freq is negative
ProbR1MinusRef = bsxfun(@minus,ProbR1Ref,probR1);

PumpR1Logic = eq(PumpR1MinusRef,zeros(size(PumpR1MinusRef)));
ProbR1Logic = eq(ProbR1MinusRef,zeros(size(ProbR1MinusRef)));

ProbAndResponseR1 = bsxfun(@times,ProbR1Logic,ResponseR1);
PumpR1Logic=+PumpR1Logic; % get ride of the logic tag

SpecAccuR1 = PumpR1Logic*ProbAndResponseR1;
% NormSpecAccuR1 = SpecAccuR1./max(abs(SpecAccuR1(:)));

%% R2
[~,PumpR2Ref]= ndgrid(ones(NumResR2,1),FreqRange);
[~,ProbR2Ref]= ndgrid(ones(NumResR2,1),FreqRange);

PumpR2MinusRef = bsxfun(@plus,PumpR2Ref,pumpR2)'; % since the pump freq is negative
ProbR2MinusRef = bsxfun(@minus,ProbR2Ref,probR2);

PumpR2Logic = eq(PumpR2MinusRef,zeros(size(PumpR2MinusRef)));
ProbR2Logic = eq(ProbR2MinusRef,zeros(size(ProbR2MinusRef)));

ProbAndResponseR2 = bsxfun(@times,ProbR2Logic,ResponseR2);
PumpR2Logic=+PumpR2Logic; % get ride of the logic tag

SpecAccuR2 = PumpR2Logic*ProbAndResponseR2;
% NormSpecAccuR2 = SpecAccuR2./max(abs(SpecAccuR2(:)));

%% R3
[~,PumpR3Ref]= ndgrid(ones(NumResR3,1),FreqRange);
[~,ProbR3Ref]= ndgrid(ones(NumResR3,1),FreqRange);

PumpR3MinusRef = bsxfun(@plus,PumpR3Ref,pumpR3)'; % since the pump freq is negative
ProbR3MinusRef = bsxfun(@minus,ProbR3Ref,probR3);

PumpR3Logic = eq(PumpR3MinusRef,zeros(size(PumpR3MinusRef)));
ProbR3Logic = eq(ProbR3MinusRef,zeros(size(ProbR3MinusRef)));

ProbAndResponseR3 = bsxfun(@times,ProbR3Logic,ResponseR3);
PumpR3Logic=+PumpR3Logic; % get ride of the logic tag

SpecAccuR3 = PumpR3Logic*ProbAndResponseR3;
% NormSpecAccuR3 = SpecAccuR3./max(abs(SpecAccuR3(:)));

%% NR1
[~,PumpNR1Ref]= ndgrid(ones(NumResNR1,1),FreqRange);
[~,ProbNR1Ref]= ndgrid(ones(NumResNR1,1),FreqRange);

PumpNR1MinusRef = bsxfun(@minus,PumpNR1Ref,pumpNR1)'; % since the pump freq is positive
ProbNR1MinusRef = bsxfun(@minus,ProbNR1Ref,probNR1);

PumpNR1Logic = eq(PumpNR1MinusRef,zeros(size(PumpNR1MinusRef)));
ProbNR1Logic = eq(ProbNR1MinusRef,zeros(size(ProbNR1MinusRef)));

ProbAndResponseNR1 = bsxfun(@times,ProbNR1Logic,ResponseNR1);
PumpNR1Logic=+PumpNR1Logic; % get ride of the logic tag

SpecAccuNR1 = PumpNR1Logic*ProbAndResponseNR1;
% NormSpecAccuNR1 = SpecAccuNR1./max(abs(SpecAccuNR1(:)));

%% NR2
[~,PumpNR2Ref]= ndgrid(ones(NumResNR2,1),FreqRange);
[~,ProbNR2Ref]= ndgrid(ones(NumResNR2,1),FreqRange);

PumpNR2MinusRef = bsxfun(@minus,PumpNR2Ref,pumpNR2)'; % since the pump freq is positive
ProbNR2MinusRef = bsxfun(@minus,ProbNR2Ref,probNR2);

PumpNR2Logic = eq(PumpNR2MinusRef,zeros(size(PumpNR2MinusRef)));
ProbNR2Logic = eq(ProbNR2MinusRef,zeros(size(ProbNR2MinusRef)));

ProbAndResponseNR2 = bsxfun(@times,ProbNR2Logic,ResponseNR2);
PumpNR2Logic=+PumpNR2Logic; % get ride of the logic tag

SpecAccuNR2 = PumpNR2Logic*ProbAndResponseNR2;
% NormSpecAccuNR2 = SpecAccuNR2./max(abs(SpecAccuNR2(:)));

%% NR3
[~,PumpNR3Ref]= ndgrid(ones(NumResNR3,1),FreqRange);
[~,ProbNR3Ref]= ndgrid(ones(NumResNR3,1),FreqRange);

PumpNR3MinusRef = bsxfun(@minus,PumpNR3Ref,pumpNR3)'; % since the pump freq is positive
ProbNR3MinusRef = bsxfun(@minus,ProbNR3Ref,probNR3);

PumpNR3Logic = eq(PumpNR3MinusRef,zeros(size(PumpNR3MinusRef)));
ProbNR3Logic = eq(ProbNR3MinusRef,zeros(size(ProbNR3MinusRef)));

ProbAndResponseNR3 = bsxfun(@times,ProbNR3Logic,ResponseNR3);
PumpNR3Logic=+PumpNR3Logic; % get ride of the logic tag

SpecAccuNR3 = PumpNR3Logic*ProbAndResponseNR3;
% NormSpecAccuNR3 = SpecAccuNR3./max(abs(SpecAccuNR3(:)));

%% Sum & output
Rephasing = SpecAccuR1 + SpecAccuR2 - SpecAccuR3;
NonRephasing = SpecAccuNR1 + SpecAccuNR2 - SpecAccuNR3;

Output.Rephasing=Rephasing;
Output.NonRephasing=NonRephasing;

Output.SpecAccuR1 =SpecAccuR1;
Output.SpecAccuR2 =SpecAccuR2;
Output.SpecAccuR3 =SpecAccuR3;
Output.SpecAccuNR1=SpecAccuNR1;
Output.SpecAccuNR2=SpecAccuNR2;
Output.SpecAccuNR3=SpecAccuNR3;

