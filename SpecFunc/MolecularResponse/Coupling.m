function [Beta,CouplingList] = Coupling(StrucInfo,CoupleType)
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

%% List for COSMOSS GUI and matching popmenu value-name pair
CouplingList = {'TDC',...
                'NN_Mix_TDC',...
                'NN_Mix_TDC_Betasheet',...
                'Cho_PB',...
                'Cho_APB',...
                'TDC+Cho_APB',...
                };

%% select models
if isstruct(StrucInfo)
    switch CoupleType
        case 'TDC'
            Beta = Coupling_TDC(StrucInfo);

        case 'NN_Mix_TDC'
            N_Modes = StrucInfo.Num_Modes;
            Beta    = Coupling_TDC(StrucInfo);
            Beta_NN = bsxfun(@times,ones(N_Modes-1,1),StrucInfo.Beta_NN);
            
            Beta(logical(diag(ones(N_Modes-1,1), 1))) = Beta_NN;
            Beta(logical(diag(ones(N_Modes-1,1),-1))) = Beta_NN;
            
        case 'NN_Mix_TDC_Betasheet'
            % this is for betasheet only
            % othe structure do not have N_residue and N_Strand
            Beta      = Coupling_TDC(StrucInfo);
            Beta_NN   = StrucInfo.Beta_NN;
            N_Residue = StrucInfo.N_Residue;
            N_Strand  = StrucInfo.N_Strand;

            N_Mode_per_Starnd = N_Residue -1;

            % Subsituting nearest neighbor Beta_NN
            for i1 = 1:N_Strand
                for j1 = 1:N_Mode_per_Starnd-1

                    Ind12 = (i1-1)*N_Mode_per_Starnd+j1;

                    Beta(Ind12    ,Ind12 + 1) = Beta_NN;
                    Beta(Ind12 + 1,Ind12    ) = Beta_NN;
                end
            end

        case 'Cho_PB'
            % this is for betasheet only
            Beta = Coupling_Cho_PB(StrucInfo);

        case 'Cho_APB'
            % this is for betasheet only
            Beta = Coupling_Cho_APB(StrucInfo);

        case 'TDC+Cho_APB'
            % this is for betasheet only
            Beta_TDC_all = Coupling_TDC(StrucInfo);
            Num_Modes1 = StrucInfo.StrucData1.Num_Modes;
            Beta_APB = Coupling_Cho_APB(StrucInfo.StrucData2);

            Beta_Mix = Beta_TDC_all;
            Beta_APB_Ind = Num_Modes1+1:StrucInfo.Num_Modes;
            Beta_Mix(Beta_APB_Ind,Beta_APB_Ind) = Beta_APB;

            Beta = Beta_Mix;
        otherwise
            disp(['Unsupported coupling model:' CoupleType])  
    end
else
    Beta = [];
end