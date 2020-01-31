function Output = MuAlphaGen(SData,H,ExMode,varargin)
%% MuAlphaGen 
% 
% This Script generate mu and alpha matrix in local mode basis. According
% to Jenny's Mathematica code, to generate alpha(Raman tensor) in local
% mode basis, the transition-slection rule is the same as mu is. 
% 
% Since the Raman tensor is symmetric for normal molecule, we can reduce
% the unique elements from 9 to 6 in GetAmideI.m. As a result, the input of
% alpha_in_Simulation_frame will be a N_modes * 6 matrix.
% 
% Todo: 1. fix [Improve] part
% 
% Copyright Jia-Jung Ho, 2013-2020

%% Inputs parser
INPUT = inputParser;

% Default values
defaultMode   = 'Mu';

expectedModes = {'Mu','Alpha'};

% Add optional inputs to inputparser object

addParamValue(INPUT,'Mode',defaultMode,...
                 @(x) any(validatestring(x,expectedModes)));
           
parse(INPUT,varargin{:});

% Reassign Variable names
Mode = INPUT.Results.Mode;

switch Mode   
    case 'Mu'
       Trans_Moment = SData.Scaled_LocMu; % size [N x 3]
       
    case 'Alpha'  
       Trans_Moment = SData.Scaled_LocAlpha; % note: RamanV = [N x 9], index: [xx xy xz yx yy yz zx zy zz]
end

M = size(Trans_Moment,2);

%% reassign variable names
N         = SData.Nmodes;
LocFreq   = SData.LocFreq;
LocAnharm = SData.LocAnharm;
Beta      = H.Beta;
OneExH    = H.H;


%% Generate Two Exciton block diag part of full Hamiltonianif needed (TwoExOvertoneH & TwoExCombinationH)
if strcmp(ExMode,'TwoEx')
    StatesNum = (N+2)*(N+1)/2; 

    %- Pre-allocate the total matrix size
    TwoExPart = zeros(StatesNum-1-N,StatesNum-1-N);

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

    % Bolck-Diagnonal Part of of TwoExCombinationH
    for m=1:N-1
        TwoExPart(TEDIndexBegin(m):TEDIndexEnd(m),TEDIndexBegin(m):TEDIndexEnd(m))...
            = bsxfun(@times,eye(N-m),LocFreq(m+1:end)+LocFreq(m)) + Beta(m+1:end,m+1:end);
    end

    %% Two Exciton Hamiltonian Cross part between TwoExOvertoneH and TwoExCombinationH
    for n1=1:N-1
        TempOffDiagMatrix = zeros(N,N-n1);
        TempOffDiagMatrix(n1,:) = Beta(n1,n1+1:N).*sqrt(2);
        TempOffDiagMatrix(n1+1:N,:) = bsxfun(@times,eye(N-n1),Beta(n1,n1+1:N).*sqrt(2));

        TwoExPart(1:N,TEDIndexBegin(n1):TEDIndexEnd(n1))...
            = TempOffDiagMatrix;
    end

    for n2=1:N-1
        TempOffDiagMatrix = zeros(N-n2,N);
        TempOffDiagMatrix(:,n2) = Beta(n2+1:N,n2).*sqrt(2);
        TempOffDiagMatrix(:,n2+1:N) = bsxfun(@times,eye(N-n2),Beta(n2,n2+1:N).*sqrt(2));

        TwoExPart(TEDIndexBegin(n2):TEDIndexEnd(n2),1:N)...
            = TempOffDiagMatrix;
    end

end

%% Diagonalize the full hamiltonian

if strcmp(ExMode,'TwoEx')
    FullH = blkdiag(OneExH,TwoExPart);
else
    FullH = OneExH;
end 

% note: the eiganvector V_Full(:,i) has been already normalized.
[V_Full,D_Full] = eig(FullH);
Ex_Freq = diag(D_Full);

% sort eiganvalue form small to big and reorder the eiganvectors
[Sort_Ex_Freq,Indx] = sort(Ex_Freq);
 Sort_Ex_V          = V_Full(:,Indx);
 
Sort_Ex_F1 = Sort_Ex_Freq(2:N+1);
Sort_Ex_V1 = Sort_Ex_V(2:N+1,2:N+1);
Sort_Ex_F2 = [];
Sort_Ex_V2 = [];

if strcmp(ExMode,'TwoEx')
    Sort_Ex_F2 = Sort_Ex_Freq(N+2:end);
    Sort_Ex_V2 = Sort_Ex_V(N+2:end,N+2:end);
end 

% export
Output.Sort_Ex_F1 = Sort_Ex_F1;
Output.Sort_Ex_V1 = Sort_Ex_V1;
Output.Sort_Ex_F2 = Sort_Ex_F2;
Output.Sort_Ex_V2 = Sort_Ex_V2;

% Zero to One exciton transition 
M_Lo_01 = Trans_Moment; % size [N x 3 or 9]

%% Two Exciton part in local mode basis
M_Lo_12 = [];

if strcmp(ExMode,'TwoEx')

    % One exciton to Overtone transition 
    Trans12_Onsite = zeros(N,N,M);

    for ii=1:N
        Trans12_Onsite(ii,ii,:) = Trans_Moment(ii,:)*sqrt(2);
    end

    % Build up the block begin/end list from Two Exciton index
    TEDIndexBegin = TEDIndexBegin-N;
    TEDIndexEnd   = TEDIndexEnd-N;

    % One exciton to Combination transition
    Trans12_Combination = zeros(N,N*(N-1)/2,M);
    for jj1=1:N-1
        Temp_Trans12_Combination          = zeros(N,N-jj1,M);
        Temp_Trans12_Combination(jj1,:,:) = Trans_Moment(jj1+1:end,:);

        % [Improve] Can remove this for loop by work out the exact indexing of 
        % the diagnoal terms.
        for kk=1:N-jj1
            Temp_Trans12_Combination(kk+jj1,kk,:) = Trans_Moment(jj1,:);
        end

        Trans12_Combination(:,TEDIndexBegin(jj1):TEDIndexEnd(jj1),:) = Temp_Trans12_Combination;
    end

    M_Lo_12 = [Trans12_Onsite,Trans12_Combination];
end

%% Change of basis and output
% 0->1, size = [N x 3 or 9] 
M_Ex_01    = Sort_Ex_V1' * M_Lo_01;

% 1->2
M_Ex_12 = [];

if strcmp(ExMode,'TwoEx')
    M_Ex_12    = zeros(N,N*(N+1)/2,M);
    
    for II = 1:M
        % size = [N x N*(N+1)/2 x 3 or 9] 
        M_Ex_12(:,:,II) = Sort_Ex_V1' * squeeze(M_Lo_12(:,:,II)) * Sort_Ex_V2; 
    end
end

% export
Output.M_Lo_01 = M_Lo_01;
Output.M_Ex_01 = M_Ex_01;
Output.M_Lo_12 = M_Lo_12;
Output.M_Ex_12 = M_Ex_12;

%% Calculate the norm or eigenvalues of the exciton modes
M_Ex_01_N = [];
M_Ex_12_N = [];
DX_01     = [];
VX_01     = [];
DX_12     = [];
VX_12     = [];

switch Mode   
    case 'Mu'
        M_Ex_01_N = sqrt(sum(M_Ex_01.^2,2)); % size = [N,1]
       
    case 'Alpha'  
        % note: RamanV = [N x 9], index: [xx xy xz yx yy yz zx zy zz]
        X_01 = reshape(M_Ex_01,[],3,3);
        [VX_01,DX_01] = eig3(X_01); 
        %M_Ex_01_N = max(abs(DX_01),[],2); % size = [N,1]
        M_Ex_01_N = sum(abs(DX_01),2); % size = [N,1], take trace
        
end

if strcmp(ExMode,'TwoEx')
    switch Mode   
        case 'Mu'
            X_12 = reshape(M_Ex_12,[],3);
            M_Ex_12_N = sqrt(sum(X_12.^2,2)); % size = [N,1]
            
        case 'Alpha'
            X_12 = reshape(M_Ex_12,[],3,3);
            [VX_12,DX_12] = eig3(X_12);
            %M_Ex_12_N = max(abs(DX_12),[],2); % size = [N,1]
            M_Ex_12_N = sum(abs(DX_12),2); % size = [N,1], take trace
            
    end
end


% Output
Output.M_Ex_01_N = M_Ex_01_N; % size = [N,1]
Output.M_Ex_12_N = M_Ex_12_N; % size = [N,1]
Output.M_Ex_01_D = DX_01; % size = [N,3]
Output.M_Ex_01_V = VX_01; % size = [N,3,3]
Output.M_Ex_12_D = DX_12; % size = [N,3]
Output.M_Ex_12_V = VX_12; % size = [N,3,3]  