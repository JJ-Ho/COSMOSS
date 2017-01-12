function AccuGrid = Bin1D(Freq,Signal,FreqRange)

    FreqMin = FreqRange(1);
    NumFreq = length(FreqRange);

    % round to 1cm-1 resolution and shift the frequency to matrix subscript
    FreqInd = round(Freq) - FreqMin + 1;
    
        
    % remove any response out of frequency range
    FreqSubMax = FreqRange(end) - FreqMin + 1;
    Freq_OutOfRange = or((FreqInd - FreqSubMax) > 0,FreqInd <=0);
    
    FreqInd(Freq_OutOfRange) = [];
    Signal(Freq_OutOfRange)  = [];
      
    AccuGrid = zeros(NumFreq,1);
    
    for i = 1:length(Signal)
        AccuGrid(FreqInd(i)) = AccuGrid(FreqInd(i)) + Signal(i);
    end



