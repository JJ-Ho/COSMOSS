function Beta = Coupling_Cho_APB(S)
%% debug
% clear all
% 
% S.N_Residue = 10;
% S.N_Strand  = 5;
% S.N_Mode_per_Starnd = S.N_Residue - 1;
% S.Num_Modes = S.N_Strand*S.N_Mode_per_Starnd;

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
StrandDirectionMatrix = ones(N_Strand,N_Mode_per_Starnd);

Residue_Array =          1:N_Mode_per_Starnd ;
Residue_Even  =        ~mod(Residue_Array,2) ;
Residue_Odd   = logical(mod(Residue_Array,2));

Strand_Array =                  1:N_Strand ;
Strand_Even  =        ~mod(Strand_Array,2) ;
Strand_Odd   = logical(mod(Strand_Array,2));

Strand1               = ones(1,N_Mode_per_Starnd);
Strand1(Residue_Even) = 1i; 
Strand2               = ones(1,N_Mode_per_Starnd).*(-1);
Strand2(Residue_Odd)  = -1i;

StrandDirectionMatrix( Strand_Odd,:) = bsxfun(@times,StrandDirectionMatrix( Strand_Odd ,:),Strand1);
StrandDirectionMatrix(Strand_Even,:) = bsxfun(@times,StrandDirectionMatrix( Strand_Even,:),Strand2);

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

Sub_Left_Direct  = StrandDirection(Sub_Left);
Sub_Right_Direct = StrandDirection(Sub_Right);

%% ap1
Pattern_ap1   = bsxfun(@times,ones(size(Sub_Left)),[-1,1]);
Match_Ind_ap1 = all(eq(Pattern_ap1,Sub_Left_Direct),2);
Sub_ap1 = Sub_Left(Match_Ind_ap1,:);
Ind_ap1 = sub2ind(size(Beta),Sub_ap1(:,1),Sub_ap1(:,2));

%% ap2
Pattern_ap2   = bsxfun(@times,ones(size(Sub_Right)),[-1i,1i]);
Match_Ind_ap2 = all(eq(Pattern_ap2,Sub_Right_Direct),2);
Sub_ap2 = Sub_Right(Match_Ind_ap2,:);
Ind_ap2 = sub2ind(size(Beta),Sub_ap2(:,1),Sub_ap2(:,2));

%% hh
Pattern_hh1   = bsxfun(@times,ones(size(Sub_Right)),[1,-1]);
Match_Ind_hh1 = all(eq(Pattern_hh1,Sub_Right_Direct),2);
Pattern_hh2   = bsxfun(@times,ones(size(Sub_Left )),[-1i,1i]);
Match_Ind_hh2 = all(eq(Pattern_hh2,Sub_Left_Direct ),2);

Sub_hh = [Sub_Right(Match_Ind_hh1,:);Sub_Left(Match_Ind_hh2,:)];
Ind_hh = sub2ind(size(Beta),Sub_hh(:,1),Sub_hh(:,2));

%% tt
Pattern_tt1   = bsxfun(@times,ones(size(Sub_Right)),[-1,1]);
Match_Ind_tt1 = all(eq(Pattern_tt1,Sub_Right_Direct),2);
Pattern_tt2   = bsxfun(@times,ones(size(Sub_Left )),[1i,-1i]);
Match_Ind_tt2 = all(eq(Pattern_tt2,Sub_Left_Direct ),2);

Sub_tt = [Sub_Right(Match_Ind_tt1,:);Sub_Left(Match_Ind_tt2,:)];
Ind_tt = sub2ind(size(Beta),Sub_tt(:,1),Sub_tt(:,2));

%% Sub Beta
Beta(Ind_ap1) = Fap1;
Beta(Ind_ap2) = Fap2;
Beta(Ind_hh)  = Fhh;
Beta(Ind_tt)  = Ftt;

%% Symmetrize
Beta = Beta + Beta';





