function Output= ExcitonH(PDB_Data,varargin)
%% TwoExcitonH
%  
% This Script generate the state transition list. I intentionally make this
% script as non-function, since there are some large matrixies, such as
% alpha_Exciton that I don't want to mvove around memory.
% 
% Todo:: 1. Should change as function input
%        2. fix [Improve] part
% 
% ------- Version log -----------------------------------------------------
%
% Ver. 3.0  140608  Unify one Ex /Two Ex H generation 
%                   (Branch out from the ExcitonH.m in 2DSFG_ssDNA)
% 
% Ver. 2.0  130901  recover the 9 element version alpha tensor for
%                   Rotational avg matrx
% 
% Ver. 1.9  130813  Disable "Generate Local States" part since it's not
%                   in use.
% 
% Ver. 1.8  130802  Fixed a bug in "Two Exciton Hamiltonian Cross part 
%                   between TwoExOnSiteH and TwoExDiffuseH". The diag(-1)
%                   terms should the same as the first row.
% 
% Ver. 1.7  130723  To accomodate the reduction of Raman tensor, the final
%                   exciton basis transformation was reduced from 9
%                   for-loops to 6 loops.
% 
% Ver. 1.6  130706  Move the FTIR ploting part to upper level =>
%                   TwoDSFG_Simulation.
% 
% Ver. 1.5  130705  Isolated the TwoDSFG_MuAlphaGen out from the original code.
% 
% Ver. 1.4  130617  Finish mu_Loc and mu_Exciton, some [Improve] can be
%                   clean up.
% 
% Ver. 1.3  130614  Finish Two Exciton Hamiltonian Cross part and
%                   correction for multipling sqrt(2) factor
% 
% Ver. 1.2  130612  Finish two exction hamiltonian and diagonalization part
% 
% Ver. 1.1  130611  Finish "OneExcitonH"
% 
% Ver. 1.0  130609  Finish "Generate Local States" part
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

%% Debug
% clear all
% PDB_Data = GetAmideI('Debug');
% ExMode = 'TwoEx';
% Coupling = 'TDC';

%% Inputs parser
INPUT = inputParser;

% Default values
defaultExMode    = 'TwoEx';
expectedExMode   = {'OneEx','TwoEx'};
defaultCoupling  = 'TDC'; 
expectedCoupling = {'TDC','Zero','NN','NN_Mix_TDC'};
defaultBeta_NN   = 0.8; % 0.8 cm-1 according to Lauren's PNAS paper (doi/10.1073/pnas.1117704109); that originate from Min Cho's paper (doi:10.1063/1.1997151)

addParamValue(INPUT,'ExMode',defaultExMode,...
                 @(x) any(validatestring(x,expectedExMode)));
addParamValue(INPUT,'Coupling',defaultCoupling,...
                 @(x) any(validatestring(x,expectedCoupling)));
addOptional(INPUT,'Beta_NN',defaultBeta_NN);
             
             
parse(INPUT,varargin{:});
             
ExMode   = INPUT.Results.ExMode;
Coupling = INPUT.Results.Coupling;
Beta_NN  = INPUT.Results.Beta_NN;

%% Initiation part
Num_Modes = PDB_Data.Num_Modes;
Center    = PDB_Data.center;
Mu        = PDB_Data.mu;
Freq      = PDB_Data.freq;
Anharm    = PDB_Data.anharm;

%% Define total State number
if strcmp(ExMode,'TwoEx')
    StatesNum = (Num_Modes+2)*(Num_Modes+1)/2; 
else
    StatesNum = (Num_Modes+1); 
end

%% Zero exciton Hamiltonain 
ZeroExH = 0;

%% One Exciton Hamiltonian
% Coupling terms in one exciton Hamiltonian

[StateIndex1,StateIndex2] = ndgrid(1:Num_Modes); % SI1 => I; SI2 => J

ModeDistVec = Center(StateIndex1(:),:) - Center(StateIndex2(:),:); % aka R
ModeDist    = sqrt(sum(abs(ModeDistVec).^2,2)); % Got to have abs, or otherwise will have imaginary part

muIdotJ = dot(Mu(StateIndex1(:),:),Mu(StateIndex2(:),:),2);
RdotmuI = dot(ModeDistVec,Mu(StateIndex1(:),:),2);
RdotmuJ = dot(ModeDistVec,Mu(StateIndex2(:),:),2);

% BetaPreFactor = 5034*((4.1058/sqrt(1600))*3.144)^2; %From Jenny's mathematica code "TransitionDipoleCoupling.m" The result is in cm-1 unit
BetaPreFactor = 84862/Freq(1)*(0.23344)^2; %From Chem. Phys. Lett. v106, #6, p613
Beta = BetaPreFactor*(muIdotJ./(ModeDist.^3) - 3*RdotmuJ.*RdotmuI./(ModeDist).^5); % Beta in Num_Modes^2 x 3 size
Beta(isnan(Beta)) = 0; % Get ride of diagnal 
Beta = reshape(Beta,Num_Modes,Num_Modes);


switch Coupling
    case 'Zero'
        Beta = zeros(size(Beta));
        
    case 'NN_Mix_TDC'
        NN_U_Index = sub2ind(size(Beta), 1:Num_Modes-1, 2:Num_Modes  );
        NN_L_Index = sub2ind(size(Beta), 2:Num_Modes  , 1:Num_Modes-1);
        Beta([NN_U_Index;NN_L_Index]) = Beta_NN;
        
    case 'NN'
        Beta = zeros(size(Beta));
        NN_U_Index = sub2ind(size(Beta), 1:Num_Modes-1, 2:Num_Modes  );
        NN_L_Index = sub2ind(size(Beta), 2:Num_Modes  , 1:Num_Modes-1);
        Beta([NN_U_Index;NN_L_Index]) = Beta_NN;
end

% The result is in cm-1 unit
OneExH = bsxfun(@times,eye(Num_Modes),Freq) + Beta;

%% Two Exciton Hamiltonian block diag part (TwoExOnSiteH & TwoExDiffuseH)
if strcmp(ExMode,'TwoEx')

    %- Pre-allocate the total matrix size
    TwoExH = zeros(StatesNum-1-Num_Modes,StatesNum-1-Num_Modes);

    %- Two Exciton Overtone Hamiltonian, TwoExOvertoneH
    TwoExH(1:Num_Modes,1:Num_Modes) = diag(2*Freq-Anharm);

    %- Two Exciton Combination Hamiltonian, TwoExCombinationH
    NumOfElementInBolcks = Num_Modes-1:-1:1;
    TEDIndexEnd   = zeros(Num_Modes-1,1);
    TEDIndexBegin = zeros(Num_Modes-1,1);

    for L=1:Num_Modes-1
        TEDIndexEnd(L) = sum(NumOfElementInBolcks(1:L)); % TED = TwoExCombinationH
        TEDIndexBegin(L) = TEDIndexEnd(L) - NumOfElementInBolcks(L) + 1;
    end

    % Shift index for accomdation of TwoExOvertoneH
    TEDIndexBegin = TEDIndexBegin+Num_Modes;
    TEDIndexEnd = TEDIndexEnd+Num_Modes;

    % Off-block-Diagonal, Lower triangular part, of  TwoExCombinationH
    for k2=1:Num_Modes-1
        for k1=k2+1:Num_Modes-1

            TempOffDiagMatrix = zeros(Num_Modes-k1,Num_Modes-k2);
            TempOffDiagMatrix(:,k1-k2) = Beta(k1+1:end,k2);
            TempOffDiagMatrix(:,k1-k2+1:end) = eye(Num_Modes-k1)*Beta(k1,k2);

            TwoExH(TEDIndexBegin(k1):TEDIndexEnd(k1),TEDIndexBegin(k2):TEDIndexEnd(k2))...
                = TempOffDiagMatrix;
        end
    end

    % Off-block-Diagonal, Upper triangular part of TwoExCombinationH
        %     
        % --- Note ----------------------------------------------------------------
        % To run this:
        % TwoExCombinationH = TwoExCombinationH + TwoExCombinationH';  
        % is slower than running another "for"
        % --- Note ----------------------------------------------------------------

    for k2=2:Num_Modes-1
        for k1=1:k2-1

            TempOffDiagMatrix = zeros(Num_Modes-k1,Num_Modes-k2);
            TempOffDiagMatrix(k2-k1,:) = Beta(k2+1:end,k1);
            TempOffDiagMatrix(k2-k1+1:end,:) = eye(Num_Modes-k2)*Beta(k2,k1);

            TwoExH(TEDIndexBegin(k1):TEDIndexEnd(k1),TEDIndexBegin(k2):TEDIndexEnd(k2))...
                = TempOffDiagMatrix;
        end
    end

    % Bolck-Diagnonal Part of of TwoExCombinationH
    for m=1:Num_Modes-1
        TwoExH(TEDIndexBegin(m):TEDIndexEnd(m),TEDIndexBegin(m):TEDIndexEnd(m))...
            = diag(Freq(m+1:end)+Freq(m)) + Beta(m+1:end,m+1:end);
    end

    %% Two Exciton Hamiltonian Cross part between TwoExOnSiteH and TwoExDiffuseH
    for n1=1:Num_Modes-1
        TempOffDiagMatrix                   = zeros(Num_Modes,Num_Modes-n1);
        TempOffDiagMatrix(n1,:)             = Beta(n1,n1+1:Num_Modes).*sqrt(2);
        TempOffDiagMatrix(n1+1:Num_Modes,:) = bsxfun(@times,eye(Num_Modes-n1),Beta(n1,n1+1:Num_Modes).*sqrt(2));

        TwoExH(1:Num_Modes,TEDIndexBegin(n1):TEDIndexEnd(n1))...
            = TempOffDiagMatrix;
    end

    for n2=1:Num_Modes-1
        TempOffDiagMatrix                   = zeros(Num_Modes-n2,Num_Modes);
        TempOffDiagMatrix(:,n2)             = Beta(n2+1:Num_Modes,n2).*sqrt(2);
        TempOffDiagMatrix(:,n2+1:Num_Modes) = bsxfun(@times,eye(Num_Modes-n2),Beta(n2,n2+1:Num_Modes).*sqrt(2));

        TwoExH(TEDIndexBegin(n2):TEDIndexEnd(n2),1:Num_Modes)...
            = TempOffDiagMatrix;
    end
end 

%% Construct Full H

if strcmp(ExMode,'TwoEx')

    H = blkdiag(ZeroExH,OneExH,TwoExH);
else 
    H = blkdiag(ZeroExH,OneExH);
end 

%% Diagonalize The full hamiltonian
% note: the eiganvector V_Full(:,i) has been already normalized.
[V_Full,D_Full] = eig(H);
Ex_Freq = diag(D_Full);

% sort eiganvalue for small to big and reorder the eiganvectors
[Sort_Ex_Freq,Indx] = sort(Ex_Freq);
Sort_Ex_V = V_Full(:,Indx);


%% Output Variables

Output.ExMode        = ExMode;
Output.StatesNum     = StatesNum;
Output.Sort_Ex_Freq  = Sort_Ex_Freq;
Output.Sort_Ex_V     = Sort_Ex_V;
Output.Beta          = Beta;
Output.H             = H;
Output.OneExH        = OneExH;


if strcmp(ExMode,'TwoEx')
    Output.TwoExH        = TwoExH;
    Output.TEDIndexBegin = TEDIndexBegin;
    Output.TEDIndexEnd   = TEDIndexEnd;
end
 
