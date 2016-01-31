function Beta = Coupling_Cho_APB(S)
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
F12  =   1.2;
F13  =   1.7;
Fab  = -10.6;
Fac  =  -1.7;
Fap1 =  -6.9;
Fap2 =  -0.7;
Fhh  =   3.2;
Ftt  =   3.9;

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
Strand1 = ones(N_Mode_per_Starnd,1);
Strand1(~mod(1:N_Mode_per_Starnd,2)) = 1i;

StrandDirectionMatrix = zeros(N_Strand,N_Mode_per_Starnd);
StrandDirectionMatrix(1,:) = Strand1;

for k = 2:N_Strand
    StrandDirectionMatrix(k,:) = fliplr(StrandDirectionMatrix(k-1,:)*-1);
end
StrandDirection = reshape(permute(StrandDirectionMatrix,[2,1]),[],1);


I1_Down  = 1:(Num_Modes-N_Mode_per_Starnd);
L_I1     = length(I1_Down);

I2_Left  = I1_Down + N_Mode_per_Starnd -1;
I2_Right = I1_Down + N_Mode_per_Starnd +1;
 
Sub_Left  = [I1_Down;I2_Left ]';
Sub_Right = [I1_Down;I2_Right]';

% remove head or tail connecting pair
Sub_Left (I_Residue_head(I1_Down),:) = [];
Sub_Right(I_Residue_tail(I1_Down),:) = [];

P_L = StrandDirection(Sub_Left(:,1)) .*StrandDirection(Sub_Left(:,2));
P_R = StrandDirection(Sub_Right(:,1)).*StrandDirection(Sub_Right(:,2));

V_tt  = nan(L_I1,1);
V_hh  = nan(L_I1,1);
V_ap1 = nan(L_I1,1);
V_ap2 = nan(L_I1,1);

V_tt (eq(P_L, 1)) = Ftt ;
V_ap1(eq(P_L,-1)) = Fap1;
V_ap2(eq(P_R, 1)) = Fap2;
V_hh (eq(P_R,-1)) = Fhh;

V_tt (isnan(V_tt )) = [];
V_ap1(isnan(V_ap1)) = [];
V_ap2(isnan(V_ap2)) = [];
V_hh (isnan(V_hh )) = [];

Ind_tt  = sub2ind(size(Beta),Sub_Left( eq(P_L, 1),2) ,Sub_Left( eq(P_L, 1),1));
Ind_ap1 = sub2ind(size(Beta),Sub_Left( eq(P_L,-1),2) ,Sub_Left( eq(P_L,-1),1));
Ind_ap2 = sub2ind(size(Beta),Sub_Right(eq(P_L, 1),2) ,Sub_Right(eq(P_L, 1),1));
Ind_hh  = sub2ind(size(Beta),Sub_Right(eq(P_L,-1),2) ,Sub_Right(eq(P_L,-1),1));

Beta(Ind_tt)  = V_tt;
Beta(Ind_ap1) = V_ap1;
Beta(Ind_ap2) = V_ap2;
Beta(Ind_hh)  = V_hh;

%% Symmetrize
Beta = Beta + Beta';





