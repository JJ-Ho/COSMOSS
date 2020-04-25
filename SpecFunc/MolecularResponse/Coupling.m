function [Beta,CouplingList] = Coupling(obj_SD,CoupleType,Beta_NN)
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
%% debug
% Structure = Data_COSMOSS.Structure;
% CoupleType = 'Jasen_TDC';

%% List for COSMOSS GUI and matching popmenu value-name pair
CouplingList = {'TDC',...
                'NN_Mix_TDC',...
                'Jansen_TDC',...
                'Zero',...
                'Constant',...
                };

%% select models

switch CoupleType
    case 'TDC'
        Beta = Coupling_TDC(obj_SD);
        
    case 'NN_Mix_TDC'
        [Beta,DistM] = Coupling_TDC(obj_SD);

        DistCutOff = 4; %Need to export this as a GUI option
        Beta(DistM < DistCutOff) = Beta_NN;

    case 'Jansen_TDC'
        Beta = Coupling_TDC(obj_SD);
        [Beta_Jansen_NN,~] = Coupling_Jansen(obj_SD);
        
        Sub_Ind = ~eq(Beta_Jansen_NN,0);
        Beta(Sub_Ind) = Beta_Jansen_NN(Sub_Ind);
        
    case 'NN_Mix_TDC_Betasheet'
        % this is for betasheet only
        % othe structure do not have N_residue and N_Strand
        Beta      = Coupling_TDC(obj_SD);
        N_Residue = obj_SD.N_Residue;
        N_Strand  = obj_SD.N_Strand;

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
        Beta = Coupling_Cho_PB(obj_SD);

    case 'Cho_APB'
        % this is for betasheet only
        Beta = Coupling_Cho_APB(obj_SD);

    case 'TDC+Cho_APB'
        % this is for betasheet only
        Beta_TDC_all = Coupling_TDC(obj_SD);
        Nmodes1 = obj_SD.StrucData1.Nmodes;
        Beta_APB = Coupling_Cho_APB(obj_SD.StrucData2);

        Beta_Mix = Beta_TDC_all;
        Beta_APB_Ind = Nmodes1+1:obj_SD.Nmodes;
        Beta_Mix(Beta_APB_Ind,Beta_APB_Ind) = Beta_APB;

        Beta = Beta_Mix;

    case 'Zero'
        Beta = zeros(obj_SD.Nmodes);
        
    case 'Constant'
        %Beta = ones(obj_SD.Nmodes).* Beta_NN;
         Beta = (ones(obj_SD.Nmodes)-diag(ones(obj_SD.Nmodes,1))).* Beta_NN;  %Add Constant coupling term to all cross terms.
    case 'TDC_PBC'
        Beta = Coupling_TDC_PBC(obj_SD);
        
    case 'Symbolic'
        Beta = sym('B%d_%d',obj_SD.Nmodes);
        Beta(logical(eye(obj_SD.Nmodes))) = zeros(obj_SD.Nmodes,1);
        obj_SD.Beta = Beta;
        
    case 'List'
        % do nothing, this is for generating the coupling list 
        Beta = [];
        
    otherwise
        disp(['Unsupported coupling model:' CoupleType])
        Beta = [];
end
