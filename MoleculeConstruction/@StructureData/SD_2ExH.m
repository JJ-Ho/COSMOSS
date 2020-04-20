function Output = SD_2ExH(obj_SD)
%% Reassigne variable names
N            = obj_SD.Nmodes;
LocFreq      = obj_SD.LocFreq;
LocAnharm    = obj_SD.LocAnharm;
Beta         = obj_SD.Beta;

%% Check if generated in Symbolic mode
symMode = false;
if isobject(LocFreq)
    LocAnharm = sym('A%d',[N,1]);
    symMode = true;
end

%% Generate Two Exciton block diag part of full Hamiltonianif needed (TwoExOvertoneH & TwoExCombinationH)
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

%% Output
Output.TEDIndexBegin = TEDIndexBegin;
Output.TEDIndexEnd   = TEDIndexEnd;
Output.TwoExPart     = TwoExPart;



