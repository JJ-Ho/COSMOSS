function Output = MuAlphaGen(PDB_Data,H,varargin)
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
     
       Trans_Moment = PDB_Data.mu; % size [N x 3]
       
    case 'Alpha'
     
       Trans_Moment = PDB_Data.alpha; % note: RamanV = [N x 9], index: [xx xy xz yx yy yz zx zy zz]
end

Size_Trans   = size(Trans_Moment,2);

%% reassign variable names

Num_Modes     = PDB_Data.Num_Modes;
ExMode        = H.ExMode;
StatesNum     = H.StatesNum;
Sort_Ex_V     = H.Sort_Ex_V;

if strcmp(ExMode,'TwoEx')
    TEDIndexBegin = H.TEDIndexBegin;
    TEDIndexEnd   = H.TEDIndexEnd;
end


Trans_Loc = zeros(StatesNum,StatesNum,Size_Trans);

% Zero to One exciton transition
Trans_Loc(1,2:Num_Modes+1,:) = Trans_Moment(:,:);
Trans_Loc(2:Num_Modes+1,1,:) = Trans_Moment(:,:);

%% Two Exciton part
if strcmp(ExMode,'TwoEx')
    % One exciton to Overtone transition 
    Temp_Trans12_Onsite = zeros(Num_Modes,Num_Modes,Size_Trans);

    for ii=1:Num_Modes
        Temp_Trans12_Onsite(ii,ii,:) = Trans_Moment(ii,:)*sqrt(2);
    end

    Trans_Loc(2:Num_Modes+1,Num_Modes+2:2*Num_Modes+1,:) = Temp_Trans12_Onsite;
    Trans_Loc(Num_Modes+2:2*Num_Modes+1,2:Num_Modes+1,:) = Temp_Trans12_Onsite;

    % Build up the block begin/end list from Two Exciton index
    TransBlockIndex1 = TEDIndexBegin + Num_Modes + 1;
    TransBlockIndex2 = TEDIndexEnd   + Num_Modes + 1;

    % One exciton to Combination transition 
    for jj1=1:Num_Modes-1
        Temp_Trans12_Combination          = zeros(Num_Modes,Num_Modes-jj1,Size_Trans);
        Temp_Trans12_Combination(jj1,:,:) = Trans_Moment(jj1+1:end,:);

        % [Improve] Can remove this for loop by work out the exact indexing of 
        % the diagnoal terms.
        for kk=1:Num_Modes-jj1
            Temp_Trans12_Combination(kk+jj1,kk,:) = Trans_Moment(jj1,:);
        end

        Trans_Loc(2:Num_Modes+1,TransBlockIndex1(jj1):TransBlockIndex2(jj1),:) = Temp_Trans12_Combination;
    end

    for jj2=1:Num_Modes-1
        Temp_Trans12_Combination          = zeros(Num_Modes-jj2,Num_Modes,Size_Trans);
        Temp_Trans12_Combination(:,jj2,:) = Trans_Moment(jj2+1:end,:);

        % [Improve] Can remove this for loop by work out the exact indexing of 
        % the diagnoal terms.
        for kk=1:Num_Modes-jj2
            Temp_Trans12_Combination(kk,kk+jj2,:) = Trans_Moment(jj2,:);
        end

        Trans_Loc(TransBlockIndex1(jj2):TransBlockIndex2(jj2),2:Num_Modes+1,:) = Temp_Trans12_Combination;
    end
end

%% Change of basis

Trans_Ex = zeros(StatesNum,StatesNum,Size_Trans);
% [Improve], maybe able to do more sophisticated
for C_mu=1:Size_Trans
    Trans_Ex(:,:,C_mu) = Sort_Ex_V'*Trans_Loc(:,:,C_mu)*Sort_Ex_V;
end

%% Output

Output.Trans_Loc = Trans_Loc;
Output.Trans_Ex  = Trans_Ex;

