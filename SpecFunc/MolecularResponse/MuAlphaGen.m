function Output = MuAlphaGen(SData,H,varargin)
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
% ------- Version log -----------------------------------------------------
%
% Ver. 3.0  141126  Instead of calculating the full matix of Transition
%                   moments, I only export 0->1 and 1->2 part
% 
% Ver. 2.0  140608  Code clean up and fix sqrt(2) error
% 
% Ver. 1.1  130723  in order to reduce unique elements of Raman tensor from 
%                   9 to 6. I generaliz this script so it can take whatever
%                   length of vectors. 
%                   Thus, "Trans_Moment" can take:
% 
%                       [mu_x, mu_y, mu_z]
%                   or  [a_xx, a_yx, a_zx, a_xy, a_yy, a_zy, a_xz, a_yz, a_zz]
%                   or  [a_xx, a_yy, a_zz, a_xy, a_yz, a_xz]
% 
%                   And remove unnecessary input "Size_Trans".
% 
% Ver. 1.0  130705  Isolated from TwoExcitonH
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013-2016

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
N      = SData.Nmodes;
ExMode = H.ExMode;

% Zero to One exciton transition
M_Lo_01 = Trans_Moment; % size [N x 3 or 9]

%% Two Exciton part
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

%% Change of basis and output

% 0->1, size = [N x 3 or 9] 
Sort_Ex_V1 = H.Sort_Ex_V1;
M_Ex_01    = Sort_Ex_V1' * M_Lo_01;
% export
Output.M_Lo_01 = M_Lo_01;
Output.M_Ex_01 = M_Ex_01;

% 1->2
if strcmp(ExMode,'TwoEx')
    Sort_Ex_V2 = H.Sort_Ex_V2;
    M_Ex_12    = zeros(N,N*(N+1)/2,M);
    
    for II = 1:M
        % size = [N x N*(N+1)/2 x 3 or 9] 
        M_Ex_12(:,:,II) = Sort_Ex_V1' * squeeze(M_Lo_12(:,:,II)) * Sort_Ex_V2; 
    end
    % export
    Output.M_Lo_12 = M_Lo_12;
    Output.M_Ex_12 = M_Ex_12;
end

%% Calculate the norm or eigenvalues of the exciton modes
switch Mode   
    case 'Mu'
        M_Ex_01_N = sqrt(sum(M_Ex_01.^2,2)); % size = [N,1]
       
    case 'Alpha'  
        % note: RamanV = [N x 9], index: [xx xy xz yx yy yz zx zy zz]
        X_01 = reshape(M_Ex_01,[],3,3);
        [VX_01,DX_01] = eig3(X_01); 
        %M_Ex_01_N = max(abs(DX_01),[],2); % size = [N,1]
        M_Ex_01_N = sum(abs(DX_01),2); % size = [N,1], take trace
        
        Output.M_Ex_01_D = DX_01; % size = [N,3]
        Output.M_Ex_01_V = VX_01; % size = [N,3,3]    
        
end
Output.M_Ex_01_N = M_Ex_01_N; % size = [N,1]

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
            
            Output.M_Ex_12_D = DX_12; % size = [N,3]
            Output.M_Ex_12_V = VX_12; % size = [N,3,3]
            
    end
    Output.M_Ex_12_N = M_Ex_12_N; % size = [N,1]
end
