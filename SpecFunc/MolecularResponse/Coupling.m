function Beta = Coupling(S,CoupleType)
%% Coupling
%  
% Coupling.m is a function to generate mode coupling with a given model. 
% The model invole: 
%   TDC         Transition Dipole Coupling
%   NN_Mix_TDC  TDC mix with designate Nearest-Neighbor coupling
%   NN_Cho      The NN model in Cho's paper doi:10.1063/1.1997151

% 
% ------- Version log -----------------------------------------------------
%
% Ver. 1.0  141002  Isolate from ExcitonH. Add Cho's betasheet model
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014



%% select models

switch CoupleType
    case 'TDC'
        Beta = Coupling_TDC(S);
        
    case 'NN_Mix_TDC'
        Beta = Coupling_TDC(S);
        Beta_NN     = S.Beta_NN; 
        N_Residue = S.N_Residue;
        N_Strand  = S.N_Strand;

        % Subsituting nearest neighbor Beta_NN
        for i1 = 1:N_Strand
            for j1 = 1:N_Residue-1

                Ind12 = (i1-1)*N_Residue+j1; 

                Beta(Ind12    ,Ind12 + 1) = Beta_NN;
                Beta(Ind12 + 1,Ind12    ) = Beta_NN;
            end
        end
        
    case 'Cho_PB'
        Beta = Coupling_Cho_PB(S);
        
    case 'Cho_APB'
        Beta = Coupling_Cho_APB(S);
end      
