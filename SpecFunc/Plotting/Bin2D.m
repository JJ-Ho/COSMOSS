function Output = Bin2D(Response,FreqRange)
%% Bin2D
% 
% Given the response in the form [ Pump, Pump-Probe, Probe, Signal ], this 
% function sort the signal to the corresponding 2D frequecny grid.
% ------- Version log -----------------------------------------------------
% 
% Ver. 1.0  160327  This is a for loop version of 2D binning intend to
%                   solve the memory issue for lerger system.  
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2016

%% Frequency Grid related parameters


%% Rename signal
R1 = Response.BinR1;
R2 = Response.BinR2;
R3 = Response.BinR3;
NR1 = Response.BinNR1;
NR2 = Response.BinNR2;
NR3 = Response.BinNR3;

%% Sort the responses to 2D frequency grid
SpecAccuR1 = SortGrid(R1(:,1),R1(:,3),FreqRange,R1(:,4));
SpecAccuR2 = SortGrid(R2(:,1),R2(:,3),FreqRange,R2(:,4));
SpecAccuR3 = SortGrid(R3(:,1),R3(:,3),FreqRange,R3(:,4));

SpecAccuNR1 = SortGrid(NR1(:,1),NR1(:,3),FreqRange,NR1(:,4));
SpecAccuNR2 = SortGrid(NR2(:,1),NR2(:,3),FreqRange,NR2(:,4));
SpecAccuNR3 = SortGrid(NR3(:,1),NR3(:,3),FreqRange,NR3(:,4));

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

%% Grid sorting function
function AccuGrid = SortGrid(PumpFreq,ProbFreq,FreqRange,Signal)

    FreqMin = FreqRange(1);
    NumFreq = length(FreqRange);

    % round to 1cm-1 resolution and shift the frequency to matrix subscript
    PumpSub = abs(round(PumpFreq)) - FreqMin + 1;
    ProbSub = abs(round(ProbFreq)) - FreqMin + 1;
        
    % remove any response out of frequency range
    FreqSubMax = FreqRange(end) - FreqMin + 1;
    Pump_OutOfRange = or((PumpSub - FreqSubMax) > 0,PumpSub <=0);
    Prob_OutOfRange = or((ProbSub - FreqSubMax) > 0,ProbSub <=0);
    PumpProb_OutOfRange = or(Pump_OutOfRange,Prob_OutOfRange);
    
    PumpSub(PumpProb_OutOfRange) = [];
    ProbSub(PumpProb_OutOfRange) = [];
    Signal(PumpProb_OutOfRange) = [];
    
    % convert subscripts to indices
    ResponseList_Ind = sub2ind([NumFreq,NumFreq],PumpSub,ProbSub);
      
    AccuGrid = zeros(NumFreq);
    
    for i = 1:length(Signal)
        AccuGrid(ResponseList_Ind(i)) = AccuGrid(ResponseList_Ind(i)) + Signal(i);
    end



