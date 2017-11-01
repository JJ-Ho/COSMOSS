function Output= ExcitonH(Structure,GUI_Inputs,ExMode)

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
% Structure  = Data_COSMOSS.Structure;
% GUI_Inputs = ParseGUI_Main(Data_COSMOSS.hGUIs);
% ExMode     = 'OneEx'; % 'OneEx' or 'TwoEx'

%% Inputs parser
% Turn Output from Read GUI to cell array
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

defaultLocFreqType  = 1;
defaultCouplingType = 'TDC';
defaultSampling     = 0;
defaultP_FlucCorr   = 100;
% defaultDD_FWHM      = 0;
defaultODD_FWHM     = 0;
defaultBeta_NN      = 0.8; % 0.8 cm-1 according to Lauren's PNAS paper (doi/10.1073/pnas.1117704109); that originate from Min Cho's paper (doi:10.1063/1.1997151)

addOptional(INPUT,'LocFreqType' ,defaultLocFreqType);
addOptional(INPUT,'CouplingType',defaultCouplingType);
addOptional(INPUT,'Sampling'    ,defaultSampling);
addOptional(INPUT,'P_FlucCorr'  ,defaultP_FlucCorr);
% addOptional(INPUT,'DD_FWHM'     ,defaultDD_FWHM);
addOptional(INPUT,'ODD_FWHM'    ,defaultODD_FWHM);
addOptional(INPUT,'Beta_NN'     ,defaultBeta_NN);

parse(INPUT,GUI_Inputs_C{:});

% Reassign Variable names
LocFreqType  = INPUT.Results.LocFreqType;
CouplingType = INPUT.Results.CouplingType;
Sampling     = INPUT.Results.Sampling;
P_FlucCorr   = INPUT.Results.P_FlucCorr;
% DD_FWHM      = INPUT.Results.DD_FWHM;
ODD_FWHM     = INPUT.Results.ODD_FWHM;
Beta_NN      = INPUT.Results.Beta_NN;

%% para variables
% Reassign variable names from StrucInfo
Nmodes          = Structure.Nmodes;
LocFreq         = Structure.LocFreq;
LocAnharm       = Structure.LocAnharm;
DiagDisorder    = Structure.DiagDisorder;
% OffDiagDisorder = Structure.OffDiagDisorder; have not implement yet

if strcmp(ExMode,'TwoEx')
    StatesNum = (Nmodes+2)*(Nmodes+1)/2; 
else
    StatesNum = (Nmodes+1); 
end

% check if apply random sampling
if Sampling
    DD_std  = DiagDisorder./(2*sqrt(2*log(2)));
    ODD_std = ODD_FWHM/(2*sqrt(2*log(2)));
else
    DD_std  = 0;
    ODD_std = 0;
end

%% Diagonal disorder if any
if eq(LocFreqType,2)
    [~,dF_Jansen] = Coupling_Jansen(Structure);
    LocFreq = LocFreq + dF_Jansen;
end

P_FlucCorr = P_FlucCorr/100; % turn percentage to number within 0~1

Correlation_Dice = rand;
if Correlation_Dice < P_FlucCorr
    dF_DD = DD_std.*(randn(1,1).*ones(Nmodes,1));
else 
    dF_DD = DD_std.*randn(Nmodes,1); 
end
LocFreq = LocFreq + dF_DD;

%% Off diagonal disorder
dBeta   = ODD_std*randn(Nmodes);
dBeta   = (dBeta + dBeta')./2; % symetrize
Beta    = Coupling(Structure,CouplingType,Beta_NN); % Coupling
Beta    = Beta + dBeta;

%% Zero exciton part of full Hamiltonain 
ZeroExPart = 0;

%% One Exciton part of full Hamiltonian
% The result is in cm-1 unit
OneExPart = bsxfun(@times,eye(Nmodes),LocFreq) + Beta;

%% Two Exciton block diag part of full Hamiltonian (TwoExOvertoneH & TwoExCombinationH)
if strcmp(ExMode,'TwoEx')

    %- Pre-allocate the total matrix size
    TwoExPart = zeros(StatesNum-1-Nmodes,StatesNum-1-Nmodes);

    %- Two Exciton Overtone Hamiltonian, TwoExOvertoneH
    TwoExPart(1:Nmodes,1:Nmodes) = diag(2*LocFreq-LocAnharm);

    %- Two Exciton Combination Hamiltonian, TwoExCombinationH
    NumOfElementInBolcks = Nmodes-1:-1:1;
    TEDIndexEnd = zeros(Nmodes-1,1);
    TEDIndexBegin = zeros(Nmodes-1,1);

    for L=1:Nmodes-1
        TEDIndexEnd(L) = sum(NumOfElementInBolcks(1:L)); % TED =TwoExCombination
        TEDIndexBegin(L) = TEDIndexEnd(L) - NumOfElementInBolcks(L) + 1;
    end

    % Shift index for accomdation of TwoExOvertoneH
    TEDIndexBegin = TEDIndexBegin+Nmodes;
    TEDIndexEnd = TEDIndexEnd+Nmodes;

    % Off-block-Diagonal, Lower triangular part, of TwoExCombinationH
    for k2=1:Nmodes-1
        for k1=k2+1:Nmodes-1

            TempOffDiagMatrix = zeros(Nmodes-k1,Nmodes-k2);
            TempOffDiagMatrix(:,k1-k2) = Beta(k1+1:end,k2);
            TempOffDiagMatrix(:,k1-k2+1:end) = eye(Nmodes-k1)*Beta(k1,k2);

            TwoExPart(TEDIndexBegin(k1):TEDIndexEnd(k1),TEDIndexBegin(k2):TEDIndexEnd(k2))...
                = TempOffDiagMatrix;
        end
    end

    % Off-block-Diagonal, Upper triangular part of TwoExCombinationH
    for k2=2:Nmodes-1
        for k1=1:k2-1

            TempOffDiagMatrix = zeros(Nmodes-k1,Nmodes-k2);
            TempOffDiagMatrix(k2-k1,:) = Beta(k2+1:end,k1);
            TempOffDiagMatrix(k2-k1+1:end,:) = eye(Nmodes-k2)*Beta(k2,k1);

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
    for m=1:Nmodes-1
        TwoExPart(TEDIndexBegin(m):TEDIndexEnd(m),TEDIndexBegin(m):TEDIndexEnd(m))...
            = bsxfun(@times,eye(Nmodes-m),LocFreq(m+1:end)+LocFreq(m)) + Beta(m+1:end,m+1:end);
    end

    %% Two Exciton Hamiltonian Cross part between TwoExOvertoneH and TwoExCombinationH
    for n1=1:Nmodes-1
        TempOffDiagMatrix = zeros(Nmodes,Nmodes-n1);
        TempOffDiagMatrix(n1,:) = Beta(n1,n1+1:Nmodes).*sqrt(2);
        TempOffDiagMatrix(n1+1:Nmodes,:) = bsxfun(@times,eye(Nmodes-n1),Beta(n1,n1+1:Nmodes).*sqrt(2));

        TwoExPart(1:Nmodes,TEDIndexBegin(n1):TEDIndexEnd(n1))...
            = TempOffDiagMatrix;
    end

    for n2=1:Nmodes-1
        TempOffDiagMatrix = zeros(Nmodes-n2,Nmodes);
        TempOffDiagMatrix(:,n2) = Beta(n2+1:Nmodes,n2).*sqrt(2);
        TempOffDiagMatrix(:,n2+1:Nmodes) = bsxfun(@times,eye(Nmodes-n2),Beta(n2,n2+1:Nmodes).*sqrt(2));

        TwoExPart(TEDIndexBegin(n2):TEDIndexEnd(n2),1:Nmodes)...
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

% Extract block diagonals
%Sort_Ex_V0 = Sort_Ex_V(1,1);
%Output.Sort_Ex_V0 = Sort_Ex_V0;

Sort_Ex_F1 = Sort_Ex_Freq(2:Nmodes+1);
Sort_Ex_V1 = Sort_Ex_V(2:Nmodes+1,2:Nmodes+1);
Output.Sort_Ex_F1 = Sort_Ex_F1;
Output.Sort_Ex_V1 = Sort_Ex_V1;

if strcmp(ExMode,'TwoEx')
    Sort_Ex_F2 = Sort_Ex_Freq(Nmodes+2:end);
    Sort_Ex_V2 = Sort_Ex_V(Nmodes+2:end,Nmodes+2:end);
    Output.Sort_Ex_F2 = Sort_Ex_F2;
    Output.Sort_Ex_V2 = Sort_Ex_V2;
end 


% Sparse_TwoExH = sparse(blkdiag(ZeroExPart,OneExPart,TwoExPart));
% [V_S_Full,D_S_Full] = eigs(Sparse_TwoExH);

% [V_OneEx,D_OneEx] = eig(OneExPart);
% [V_TwoEx,D_TwoEx] = eig(TwoExPart);
% [V_TwoExOvertone,D_TwoExOvertone] = eig(TwoExOvertoneH);
% [V_TwoExCombination,D_TwoExCombination] = eig(TwoExCombinationH);

%% Output Variables
Output.ExMode        = ExMode;
Output.Nmodes        = Nmodes;
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