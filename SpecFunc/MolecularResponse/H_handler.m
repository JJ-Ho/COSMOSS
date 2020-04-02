function Output= H_handler(SData,Main_GUI_Inputs,ExMode)
% this function handle the One exciton Hamiltonian in the StrutureData with
% the following tasks:
%   1. Diagonal disorder
%   2. Off-Diagonal disorder
%   3. Generate Two-Exciton part if needed
%   4. Diagnalization

%% Inputs parser
% Turn Output from Read GUI to cell array
GUI_Inputs_C      = fieldnames(Main_GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(Main_GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = true;

defaultSampling     = 0;
defaultDD_FWHM      = 0;
defaultODD_FWHM     = 0;
defaultP_FlucCorr   = 100;

addOptional(INPUT,'Sampling'    ,defaultSampling);
addOptional(INPUT,'DD_FWHM'     ,defaultDD_FWHM);
addOptional(INPUT,'ODD_FWHM'    ,defaultODD_FWHM);
addOptional(INPUT,'P_FlucCorr'  ,defaultP_FlucCorr);

parse(INPUT,GUI_Inputs_C{:});

Sampling     = INPUT.Results.Sampling;
DD_FWHM      = INPUT.Results.DD_FWHM;
ODD_FWHM     = INPUT.Results.ODD_FWHM;
P_FlucCorr   = INPUT.Results.P_FlucCorr;

%% reassign variable names
OneExH       = SData.OneExH;
Beta         = SData.Beta;
N            = SData.Nmodes;
LocFreq      = SData.LocFreq;
LocAnharm    = SData.LocAnharm;

%% Check if generated in Symbolic mode
symMode = false;
if isobject(LocFreq)
    LocAnharm = sym('A%d',[N,1]);
    symMode = true;
end
%% check if apply random sampling
if Sampling
    DD_std  = DD_FWHM./(2*sqrt(2*log(2)));
    ODD_std = ODD_FWHM/(2*sqrt(2*log(2)));
else
    DD_std  = 0;
    ODD_std = 0;
end

%% Diagonal/Off-Diaginal disorder if any
P_FlucCorr = P_FlucCorr/100; % turn percentage to number within 0~1

% diagonal disorder
Correlation_Dice = rand;
if Correlation_Dice < P_FlucCorr
    dF_DD = DD_std.*(randn(1,1).*ones(N,1));
else 
    dF_DD = DD_std.*randn(N,1); 
end
dF_DD = [0;dF_DD]; % add the zero exciton part

% Off diagonal disorder
dBeta   = ODD_std*randn(N);
dBeta   = (dBeta + dBeta')./2; % symetrize
dBeta   = blkdiag(0,dBeta);

OneExH = OneExH + bsxfun(@times,eye(N+1),dF_DD) + dBeta;

%% Generate Two Exciton block diag part of full Hamiltonianif needed (TwoExOvertoneH & TwoExCombinationH)
TEDIndexBegin = [];
TEDIndexEnd   = [];
TwoExPart     = [];

if strcmp(ExMode,'TwoEx')
    StatesNum = (N+2)*(N+1)/2; 

    %- Pre-allocate the total matrix size
    TwoExPart = zeros(StatesNum-1-N);
    if symMode
        TwoExPart = sym(TwoExPart);
    end

    %- Two Exciton Overtone Hamiltonian, TwoExOvertoneH
    TwoExPart(1:N,1:N) = diag(2*LocFreq-LocAnharm);

    %- Two Exciton Combination Hamiltonian, TwoExCombinationH
    NumOfElementInBolcks = N-1:-1:1;
    TEDIndexEnd   = zeros(N-1,1);
    TEDIndexBegin = zeros(N-1,1);

    for L=1:N-1
        TEDIndexEnd(L)   = sum(NumOfElementInBolcks(1:L)); % TED =TwoExCombination
        TEDIndexBegin(L) = TEDIndexEnd(L) - NumOfElementInBolcks(L) + 1;
    end

    % Shift index for accomdation of TwoExOvertoneH
    TEDIndexBegin = TEDIndexBegin+N;
    TEDIndexEnd   = TEDIndexEnd+N;

    % Off-block-Diagonal, Lower triangular part, of TwoExCombinationH
    for k2=1:N-1
        for k1=k2+1:N-1

            TempOffDiagMatrix = zeros(N-k1,N-k2);
            if symMode
                TempOffDiagMatrix = sym(TempOffDiagMatrix);
            end
            
            TempOffDiagMatrix(:,k1-k2) = Beta(k1+1:end,k2);
            TempOffDiagMatrix(:,k1-k2+1:end) = eye(N-k1)*Beta(k1,k2);

            TwoExPart(TEDIndexBegin(k1):TEDIndexEnd(k1),TEDIndexBegin(k2):TEDIndexEnd(k2))...
                = TempOffDiagMatrix;
        end
    end

    % Off-block-Diagonal, Upper triangular part of TwoExCombinationH
    for k2=2:N-1
        for k1=1:k2-1

            TempOffDiagMatrix = zeros(N-k1,N-k2);
            if symMode
                TempOffDiagMatrix = sym(TempOffDiagMatrix);
            end
            
            TempOffDiagMatrix(k2-k1,:) = Beta(k2+1:end,k1);
            TempOffDiagMatrix(k2-k1+1:end,:) = eye(N-k2)*Beta(k2,k1);

            TwoExPart(TEDIndexBegin(k1):TEDIndexEnd(k1),TEDIndexBegin(k2):TEDIndexEnd(k2))...
                = TempOffDiagMatrix;
        end
    end

    % --- Note ----------------------------------------------------------------
    % To run this:
    % TwoExCombinationH = TwoExCombinationH + TwoExCombinationH';  
    % is slower than running another "for"
    % --- Note ----------------------------------------------------------------

    % Block-Diagnonal Part of of TwoExCombinationH
    for m=1:N-1
        TwoExPart(TEDIndexBegin(m):TEDIndexEnd(m),TEDIndexBegin(m):TEDIndexEnd(m))...
            =  diag(LocFreq(m+1:end)+LocFreq(m))+ Beta(m+1:end,m+1:end);
    end

    %% Two Exciton Hamiltonian Cross part between TwoExOvertoneH and TwoExCombinationH
    for n1=1:N-1
        TempOffDiagMatrix = zeros(N,N-n1);
        if symMode
            TempOffDiagMatrix = sym(TempOffDiagMatrix);
        end
        
        TempOffDiagMatrix(n1,:) = Beta(n1,n1+1:N).*sqrt(2);
        TempOffDiagMatrix(n1+1:N,:) = eye(N-n1).*Beta(n1,n1+1:N).*sqrt(2); %bsxfun(@times,eye(N-n1),Beta(n1,n1+1:N).*sqrt(2));

        TwoExPart(1:N,TEDIndexBegin(n1):TEDIndexEnd(n1))...
            = TempOffDiagMatrix;
    end

    for n2=1:N-1
        TempOffDiagMatrix = zeros(N-n2,N);
        if symMode
            TempOffDiagMatrix = sym(TempOffDiagMatrix);
        end
        
        TempOffDiagMatrix(:,n2) = Beta(n2+1:N,n2).*sqrt(2);
        TempOffDiagMatrix(:,n2+1:N) = eye(N-n2).*Beta(n2,n2+1:N).*sqrt(2);%bsxfun(@times,eye(N-n2),Beta(n2,n2+1:N).*sqrt(2));

        TwoExPart(TEDIndexBegin(n2):TEDIndexEnd(n2),1:N)...
            = TempOffDiagMatrix;
    end

end

%% Diagonalize the full hamiltonian
Sort_Ex_F1 = [];
Sort_Ex_V1 = [];
Sort_Ex_F2 = [];
Sort_Ex_V2 = [];

if strcmp(ExMode,'TwoEx')
    FullH = blkdiag(OneExH,TwoExPart);
else
    FullH = OneExH;
end 

% Diagonalization
if ~symMode
    % note: the eiganvector V_Full(:,i) has been already normalized.
    [V_Full,D_Full] = eig(FullH);
    Ex_Freq = diag(D_Full);

    % sort eiganvalue form small to big and reorder the eiganvectors
    [Sort_Ex_Freq,Indx] = sort(Ex_Freq);
     Sort_Ex_V          = V_Full(:,Indx);

    Sort_Ex_F1 = Sort_Ex_Freq(2:N+1);
    Sort_Ex_V1 = Sort_Ex_V(2:N+1,2:N+1);

    if strcmp(ExMode,'TwoEx')
        Sort_Ex_F2 = Sort_Ex_Freq(N+2:end);
        Sort_Ex_V2 = Sort_Ex_V(N+2:end,N+2:end);
    end 
end 
%% export
Output.Sort_Ex_F1 = Sort_Ex_F1;
Output.Sort_Ex_V1 = Sort_Ex_V1;
Output.Sort_Ex_F2 = Sort_Ex_F2;
Output.Sort_Ex_V2 = Sort_Ex_V2;
Output.TEDIndexBegin = TEDIndexBegin;
Output.TEDIndexEnd   = TEDIndexEnd;
Output.TwoExPart     = TwoExPart;
