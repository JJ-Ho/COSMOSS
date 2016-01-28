function Output= ExcitonH(StrucInfo,varargin)

%% TwoExcitonH
%  
% Given StrucInfo that generate by GetAmindeI, this script generate One or
% Two Exciton Hamiltonian and also export eigenvalues and eigenvectors for
% generating transition dipole and Raman tensor matrix. 
% 
% use "Mode" key word to select which H want to generate.
% expectedModes = {'OneEx','TwoEx'};
% 

% Todo: fix [Improve] part
%       Integrate Diagonal disorder part in.

% 
% ------- Version log -----------------------------------------------------
%
% Ver. 2.5  140131  Add Beta into output
% 
% Ver. 2.4  140131  Add Off-Diagnal Disorder
% 
% Ver. 2.3  140128  Modify nearest neighbor coupling to constant 0.8 cm-1
% 
% Ver. 2.2  140110  Integrate OneExH and TwoExH version by adding
%                   selection input "Mode"
% 
% Ver. 2.1  140110  Functionalize.
%                   Isolate MuAlphaGen
%                   Change Name: TwoExOnSiteH  => TwoExOvertoneH
%                                TwoExDiffuseH => TwoExCombinationH
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
% Ver. 1.5  130705  Isolated the MuAlphaGen out from the original code.
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
% varargin = {'OffDiagDisorder',0,'CouplingType','NN_Mix_TDC'};


%% Initiation part

% Reassign variable names from StrucInfo
Num_Modes = StrucInfo.Num_Modes;
Freq       = StrucInfo.freq;
Anharm     = StrucInfo.anharm;

%% Inputs parser
INPUT = inputParser;

% Default values
defaultExMode            = 'OneEx';
defaultCouplingType      = 'NN_Mix_TDC';
defaultOffDiagDisorder   = 0;
defaultBeta_NN   = 0.8; % 0.8 cm-1 according to Lauren's PNAS paper (doi/10.1073/pnas.1117704109); that originate from Min Cho's paper (doi:10.1063/1.1997151)

expectedCouplingType = {'NN_Mix_TDC','TDC','Cho_PB','Cho_APB'}; % TDC = Transition Dipole Coupling; NN = Nearest Neighbor.
expectedExModes = {'OneEx','TwoEx'};


% Add optional inputs to inputparser object
addParamValue(INPUT,'ExMode',defaultExMode,...
                 @(x) any(validatestring(x,expectedExModes)));
             
addParamValue(INPUT,'CouplingType',defaultCouplingType,...
                 @(x) any(validatestring(x,expectedCouplingType)));       
             
addOptional(INPUT,'OffDiagDisorder',defaultOffDiagDisorder);
addOptional(INPUT,'Beta_NN',defaultBeta_NN);

parse(INPUT,varargin{:});

% Reassign Variable names
ExMode          = INPUT.Results.ExMode;
CouplingType    = INPUT.Results.CouplingType;
OffDiagDisorder = INPUT.Results.OffDiagDisorder;
Beta_NN  = INPUT.Results.Beta_NN;

if strcmp(ExMode,'TwoEx')
    StatesNum = (Num_Modes+2)*(Num_Modes+1)/2; 
else
    StatesNum = (Num_Modes+1); 
end

%% Zero exciton part of full Hamiltonain 
ZeroExPart = 0;

%% One Exciton part of full Hamiltonian

StrucInfo.Beta_NN = Beta_NN; 
Beta = Coupling(StrucInfo,CouplingType);

% off diagonal disorder
OffDisorder = sqrt(OffDiagDisorder)*randn(Num_Modes);
OffDisorder = (OffDisorder + OffDisorder')./2; % symetrize
Beta = Beta + OffDisorder;

% The result is in cm-1 unit
OneExPart = bsxfun(@times,eye(Num_Modes),Freq) + Beta;

%% Two Exciton block diag part of full Hamiltonian (TwoExOvertoneH & TwoExCombinationH)
if strcmp(ExMode,'TwoEx')

    %- Pre-allocate the total matrix size
    TwoExPart = zeros(StatesNum-1-Num_Modes,StatesNum-1-Num_Modes);

    %- Two Exciton Overtone Hamiltonian, TwoExOvertoneH
    TwoExPart(1:Num_Modes,1:Num_Modes) = diag(2*Freq-Anharm);

    %- Two Exciton Combination Hamiltonian, TwoExCombinationH
    NumOfElementInBolcks = Num_Modes-1:-1:1;
    TEDIndexEnd = zeros(Num_Modes-1,1);
    TEDIndexBegin = zeros(Num_Modes-1,1);

    for L=1:Num_Modes-1
        TEDIndexEnd(L) = sum(NumOfElementInBolcks(1:L)); % TED =TwoExCombination
        TEDIndexBegin(L) = TEDIndexEnd(L) - NumOfElementInBolcks(L) + 1;
    end

    % Shift index for accomdation of TwoExOvertoneH
    TEDIndexBegin = TEDIndexBegin+Num_Modes;
    TEDIndexEnd = TEDIndexEnd+Num_Modes;

    % Off-block-Diagonal, Lower triangular part, of TwoExCombinationH
    for k2=1:Num_Modes-1
        for k1=k2+1:Num_Modes-1

            TempOffDiagMatrix = zeros(Num_Modes-k1,Num_Modes-k2);
            TempOffDiagMatrix(:,k1-k2) = Beta(k1+1:end,k2);
            TempOffDiagMatrix(:,k1-k2+1:end) = eye(Num_Modes-k1)*Beta(k1,k2);

            TwoExPart(TEDIndexBegin(k1):TEDIndexEnd(k1),TEDIndexBegin(k2):TEDIndexEnd(k2))...
                = TempOffDiagMatrix;
        end
    end

    % Off-block-Diagonal, Upper triangular part of TwoExCombinationH
    for k2=2:Num_Modes-1
        for k1=1:k2-1

            TempOffDiagMatrix = zeros(Num_Modes-k1,Num_Modes-k2);
            TempOffDiagMatrix(k2-k1,:) = Beta(k2+1:end,k1);
            TempOffDiagMatrix(k2-k1+1:end,:) = eye(Num_Modes-k2)*Beta(k2,k1);

            TwoExPart(TEDIndexBegin(k1):TEDIndexEnd(k1),TEDIndexBegin(k2):TEDIndexEnd(k2))...
                = TempOffDiagMatrix;
        end
    end

    % --- Note ----------------------------------------------------------------
    % To run this:
    % TwoExCombinationH = TwoExCombinationH + TwoExCombinationH';  
    % is slower than running another "for"
    % --- Note ----------------------------------------------------------------

    % Bolck-Diagnonal Part of of TwoExCombinationH
    for m=1:Num_Modes-1
        TwoExPart(TEDIndexBegin(m):TEDIndexEnd(m),TEDIndexBegin(m):TEDIndexEnd(m))...
            = bsxfun(@times,eye(Num_Modes-m),Freq(m+1:end)+Freq(m)) + Beta(m+1:end,m+1:end);
    end

    %% Two Exciton Hamiltonian Cross part between TwoExOvertoneH and TwoExCombinationH
    for n1=1:Num_Modes-1
        TempOffDiagMatrix = zeros(Num_Modes,Num_Modes-n1);
        TempOffDiagMatrix(n1,:) = Beta(n1,n1+1:Num_Modes).*sqrt(2);
        TempOffDiagMatrix(n1+1:Num_Modes,:) = bsxfun(@times,eye(Num_Modes-n1),Beta(n1,n1+1:Num_Modes).*sqrt(2));

        TwoExPart(1:Num_Modes,TEDIndexBegin(n1):TEDIndexEnd(n1))...
            = TempOffDiagMatrix;
    end

    for n2=1:Num_Modes-1
        TempOffDiagMatrix = zeros(Num_Modes-n2,Num_Modes);
        TempOffDiagMatrix(:,n2) = Beta(n2+1:Num_Modes,n2).*sqrt(2);
        TempOffDiagMatrix(:,n2+1:Num_Modes) = bsxfun(@times,eye(Num_Modes-n2),Beta(n2,n2+1:Num_Modes).*sqrt(2));

        TwoExPart(TEDIndexBegin(n2):TEDIndexEnd(n2),1:Num_Modes)...
            = TempOffDiagMatrix;
    end

end

%% Construct Full H
if strcmp(ExMode,'TwoEx')

    H = blkdiag(ZeroExPart,OneExPart,TwoExPart);
else 
    H = blkdiag(ZeroExPart,OneExPart);
end 

%% Diagonalize The full hamiltonian
% note: the eiganvector V_Full(:,i) has been already normalized.
[V_Full,D_Full] = eig(H);
Ex_Freq = diag(D_Full);

% sort eiganvalue form small to big and reorder the eiganvectors
[Sort_Ex_Freq,Indx] = sort(Ex_Freq);
 Sort_Ex_V          = V_Full(:,Indx);
 

% Sparse_TwoExH = sparse(blkdiag(ZeroExPart,OneExPart,TwoExPart));
% [V_S_Full,D_S_Full] = eigs(Sparse_TwoExH);

% [V_OneEx,D_OneEx] = eig(OneExPart);
% [V_TwoEx,D_TwoEx] = eig(TwoExPart);
% [V_TwoExOvertone,D_TwoExOvertone] = eig(TwoExOvertoneH);
% [V_TwoExCombination,D_TwoExCombination] = eig(TwoExCombinationH);


%% Output Variables

Output.ExMode        = ExMode;
Output.Num_Modes     = Num_Modes;
Output.StatesNum     = StatesNum;
Output.Sort_Ex_Freq  = Sort_Ex_Freq;
Output.Sort_Ex_V     = Sort_Ex_V;
Output.Beta          = Beta;
Output.H             = H;
Output.OneExPart     = OneExPart;



if strcmp(ExMode,'TwoEx')
    Output.TwoExPart     = TwoExPart;
    Output.TEDIndexBegin = TEDIndexBegin;
    Output.TEDIndexEnd   = TEDIndexEnd;
end

% Output.V_Full        = V_Full;
% Output.Indx          = Indx;
% Output.Ex_Freq       = Ex_Freq;
% Output.OneExH        = blkdiag(ZeroExPart,OneExPart);
% Output.TwoExH        = TwoExH;
