function ScoreGird = Score_Sign(SigGrid,Target)
%% Normalize each selected signal with it's max
% [S1,S2,....]./[max(S1),max(S2),...]
Norm  = max(abs(SigGrid),[],1); % max value of each column in SigGrid
Sig_N = bsxfun(@rdivide,SigGrid,Norm);  % Signal grid normalized to values with in [-1,1]

%% Pick out the targeted sign of each selected signal
Sig_N_SP = bsxfun(@times,Sig_N,sign(Target)); % for the picked sign of each Signal grid to be positive
Sig_N(Sig_N_SP<0) = 0; % force the unwanted sign of the each signal column to be zero 

%% Define score with the product of all processed signals
ScoreGird = prod(Sig_N,2);
ScoreGird = ScoreGird./max(abs(ScoreGird));
            