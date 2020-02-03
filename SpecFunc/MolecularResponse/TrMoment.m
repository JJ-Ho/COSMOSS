function Output = TrMoment(SData,ExMode,TrMode,H)
%% Transition Moment generator
% 
% This Script generate transition dipole (mu) or Raman tensor (alpha) in 
% both local mode and exciton mode basis. According to Jenny's Mathematica 
% code, to generate alpha in local mode basis, the transition-slection rule 
% is the same as mu. 
% Copyright Jia-Jung Ho, 2013-2020


%% Switch transition moment mode
switch TrMode   
    case 'Mu'
       Trans_Moment = SData.Scaled_LocMu; % size [N x 3]
       
    case 'Alpha'  
       Trans_Moment = SData.Scaled_LocAlpha; % note: RamanV = [N x 9], index: [xx xy xz yx yy yz zx zy zz]
end

M = size(Trans_Moment,2);

%% Zero to One exciton transition 
M_Lo_01 = Trans_Moment; % size [N x 3 or 9]

%% Two Exciton part in local mode basis
N = SData.Nmodes;
M_Lo_12 = [];

if strcmp(ExMode,'TwoEx')

    % One exciton to Overtone transition 
    Trans12_Onsite = zeros(N,N,M);

    for ii=1:N
        Trans12_Onsite(ii,ii,:) = Trans_Moment(ii,:)*sqrt(2);
    end

    % Build up the block begin/end list from Two Exciton index
    TEDIndexBegin = H.TEDIndexBegin-N;
    TEDIndexEnd   = H.TEDIndexEnd-N;

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

%% Conver the TrMoment from local mode to Exciton mode basis
% 0->1, size = [N x 3 or 9] 
M_Ex_01    = H.Sort_Ex_V1' * M_Lo_01;

% 1->2
M_Ex_12 = [];

if strcmp(ExMode,'TwoEx')
    M_Ex_12    = zeros(N,N*(N+1)/2,M);
    
    for II = 1:M
        % size = [N x N*(N+1)/2 x 3 or 9] 
        M_Ex_12(:,:,II) = H.Sort_Ex_V1' * squeeze(M_Lo_12(:,:,II)) * H.Sort_Ex_V2; 
    end
end

%% Calculate the norm or eigenvalues of the exciton modes
% This is for evaluating how large is such transition
M_Ex_01_N = [];
M_Ex_12_N = [];
DX_01     = [];
VX_01     = [];
DX_12     = [];
VX_12     = [];

switch TrMode   
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
    switch TrMode   
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

%% Output
Output.M_Lo_01 = M_Lo_01;
Output.M_Ex_01 = M_Ex_01;
Output.M_Lo_12 = M_Lo_12;
Output.M_Ex_12 = M_Ex_12;

% For Feynman Path Gen cutoff
Output.M_Ex_01_N = M_Ex_01_N; % size = [N,1]
Output.M_Ex_12_N = M_Ex_12_N; % size = [N,1]
Output.M_Ex_01_D = DX_01;     % size = [N,3]
Output.M_Ex_01_V = VX_01;     % size = [N,3,3]
Output.M_Ex_12_D = DX_12;     % size = [N,3]
Output.M_Ex_12_V = VX_12;     % size = [N,3,3]  