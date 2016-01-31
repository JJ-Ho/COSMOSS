function Beta = Coupling_Cho_PB(S)
%% debug
% clear all
% 
% S.N_Residue = 10;
% S.N_Strand  = 5;
% S.Num_Modes = S.N_Strand*S.N_Residue;

%% Reassign variable names from StrucInfo
Num_Modes         = S.Num_Modes;
N_Mode_per_Starnd = S.N_Mode_per_Starnd;
N_Strand          = S.N_Strand;

% Define Parallel Betasheet parameters
F12  =  0.8;
F13  =  1.7;
Fab  = -8.8;
Fac  = -1.5;
Fpb1 = -1.1;
Fpb2 =  1.5;
Fpb3 = -1.0;
Fpb4 =  1.5;

Beta = zeros(Num_Modes);

% define mode index
I_Mode = (1: Num_Modes)';
[I_Residue,~] = ind2sub([N_Mode_per_Starnd,N_Strand],I_Mode);

I_Residue_head = eq(I_Residue,1);
I_Residue_tail = eq(I_Residue,N_Mode_per_Starnd);

%% F12
Sub12 = [I_Mode,I_Mode+1];

Remove_I_12 = I_Residue_tail;
Sub12(Remove_I_12,:) = [];
Ind12 = sub2ind(size(Beta),Sub12(:,2),Sub12(:,1));
Beta(Ind12) = F12;

%% F13
Sub13 = [I_Mode,I_Mode+2];

Remove_I_13 = or(I_Residue_tail, circshift(I_Residue_tail,[-1,0]) );
Sub13(Remove_I_13,:) = [];
Ind13 = sub2ind(size(Beta),Sub13(:,2),Sub13(:,1));
Beta(Ind13) = F13;

%% Fab
Sub_ab = [I_Mode(1:Num_Modes-N_Mode_per_Starnd),I_Mode(N_Mode_per_Starnd+1:Num_Modes)];
Ind_ab = sub2ind(size(Beta),Sub_ab(:,2),Sub_ab(:,1));
Beta(Ind_ab) = Fab;

%% Fac
Sub_ac = [I_Mode(1:Num_Modes-2*N_Mode_per_Starnd),I_Mode(2*N_Mode_per_Starnd+1:Num_Modes)];
Ind_ac = sub2ind(size(Beta),Sub_ac(:,2),Sub_ac(:,1));
Beta(Ind_ac) = Fac;


%% Fpbx
I1_Down  = 1:(Num_Modes-N_Mode_per_Starnd);
L_I1     = length(I1_Down);

I2_Left  = I1_Down + N_Mode_per_Starnd -1;
I2_Right = I1_Down + N_Mode_per_Starnd +1;
 
Sub_Left  = [I1_Down;I2_Left ]';
Sub_Right = [I1_Down;I2_Right]';

Sub_pb1 = nan(L_I1,2);
Sub_pb2 = nan(L_I1,2);
Sub_pb3 = nan(L_I1,2);
Sub_pb4 = nan(L_I1,2);

for k = 1:length(I1_Down)
    if mod(I_Residue(k),2) 
        % case odd, left = Fpb4
        Sub_pb4(k,:) = Sub_Left(k,:);
        
        % case odd, right = Fpb2
        Sub_pb2(k,:) = Sub_Right(k,:);
    else
        % case even, left = Fpb1
        Sub_pb1(k,:) = Sub_Left(k,:);
        
        % case even, right = Fpb3
        Sub_pb3(k,:) = Sub_Right(k,:);
    end
end

% Count heads and tails
N_heads = sum(I_Residue_head(I1_Down));
N_tails = sum(I_Residue_tail(I1_Down));

% case left
Sub_pb1(I_Residue_head(I1_Down),:) = nan(N_heads,2);
Sub_pb4(I_Residue_head(I1_Down),:) = nan(N_heads,2);
% case right
Sub_pb2(I_Residue_tail(I1_Down),:) = nan(N_tails,2);
Sub_pb3(I_Residue_tail(I1_Down),:) = nan(N_tails,2);

% clean up NaNs
Sub_pb1(isnan(Sub_pb1(:,1)),:) = [];
Sub_pb2(isnan(Sub_pb2(:,1)),:) = [];
Sub_pb3(isnan(Sub_pb3(:,1)),:) = [];
Sub_pb4(isnan(Sub_pb4(:,1)),:) = [];

% turn sub to ind
Ind_pb1 = sub2ind(size(Beta),Sub_pb1(:,2),Sub_pb1(:,1));
Ind_pb2 = sub2ind(size(Beta),Sub_pb2(:,2),Sub_pb2(:,1));
Ind_pb3 = sub2ind(size(Beta),Sub_pb3(:,2),Sub_pb3(:,1));
Ind_pb4 = sub2ind(size(Beta),Sub_pb4(:,2),Sub_pb4(:,1));


Beta(Ind_pb1) = Fpb1;
Beta(Ind_pb2) = Fpb2;
Beta(Ind_pb3) = Fpb3;
Beta(Ind_pb4) = Fpb4;

%% Symmetrize
Beta = Beta + Beta';





