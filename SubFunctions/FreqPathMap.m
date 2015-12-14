 function R_pumpS_probe = FreqPathMap(Response,PumpMin,PumpMax,ProbeMin,ProbeMax)
% The FreqPathMap function can narror down the corresponding pathways that
% contribute to a frequency window defined by PumpRange X ProbRange

%% debug
% PumpMin = 0;
% PumpMax = 0.5;
% ProbeMin = 0;
% ProbeMax = 0.5;
% 
% Response = rand(100,5);

%% function begin
% prep
Num_Path = size(Response,1);
Ind = (1:Num_Path)';
Freq = Response(:,1:3);
Resp = Response(:,4:end);
R = [Ind abs(Freq) Resp]; % recombine the response function with path index and all positive pump/probe frequency

% extract paths in pump window
I_pump = R(:,2) > PumpMin & R(:,2) < PumpMax;
R_pump = R(I_pump,:);


% extract paths in probe window
I_probe = R_pump(:,4) > ProbeMin & R_pump(:,4) < ProbeMax;
R_pump_probe = R_pump(I_probe,:);

% sort pump frequemcy
[~,I_S] = sort(R_pump_probe(:,2));
R_pumpS_probe = R_pump_probe(I_S,:);