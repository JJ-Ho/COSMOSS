function [Beta,CouplingList] = Coupling(SData,CoupleType,Beta_NN)
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
                'Zero',...
                };

%% select models

switch CoupleType
    case 'TDC'
        Beta = Coupling_TDC(SData);

    case 'NN_Mix_TDC'
        [Beta,DistM] = Coupling_TDC(SData);

        DistCutOff = 4;
        Beta(DistM < DistCutOff) = Beta_NN;

    case 'NN_Mix_TDC_Betasheet'
        % this is for betasheet only
        % othe structure do not have N_residue and N_Strand
        Beta      = Coupling_TDC(SData);
        N_Residue = SData.N_Residue;
        N_Strand  = SData.N_Strand;

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
        Beta = Coupling_Cho_PB(SData);

    case 'Cho_APB'
        % this is for betasheet only
        Beta = Coupling_Cho_APB(SData);

    case 'TDC+Cho_APB'
        % this is for betasheet only
        Beta_TDC_all = Coupling_TDC(SData);
        Nmodes1 = SData.StrucData1.Nmodes;
        Beta_APB = Coupling_Cho_APB(SData.StrucData2);

        Beta_Mix = Beta_TDC_all;
        Beta_APB_Ind = Nmodes1+1:SData.Nmodes;
        Beta_Mix(Beta_APB_Ind,Beta_APB_Ind) = Beta_APB;

        Beta = Beta_Mix;

    case 'Zero'
        Nmodes = SData.Nmodes;
        Beta = zeros(Nmodes);
    case 'TDC_PBC'
        Beta = Coupling_TDC_PBC(SData);

    case 'None'
        % do nothing, this is for generating the coupling list 
        Beta = [];
        
    otherwise
        disp(['Unsupported coupling model:' CoupleType])
        Beta = [];
end
