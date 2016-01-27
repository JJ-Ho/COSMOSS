function Beta = Coupling_Cho_PB(S)
%% debug
% S.N_Residue = 4;
% S.N_Strand  = 3;
% S.Num_Modes = S.N_Strand*S.N_Residue;

%% Reassign variable names from StrucInfo
Num_Modes = S.Num_Modes;
N_Residue = S.N_Residue;
N_Strand  = S.N_Strand;

Beta = zeros(Num_Modes);

% Define Parallel Betasheet parameters
F12  =  0.8;
F13  =  1.7;
Fab  = -8.8;
Fac  = -1.5;
Fpb1 = -1.1;
Fpb2 =  1.5;
Fpb3 = -1.0;
Fpb4 =  1.5;

% N_pair = N_Strand*(N_Residue-1);
N_pair = Num_Modes*(Num_Modes-1)/2;

%% F12
Ind12 = nan(N_pair,2);

for i1 = 1:N_Strand
    for j1 = 1:N_Residue-1
        
        N_loop = (i1-1)*N_Residue+j1;
        Ind12(N_loop,:) = [(i1-1)*N_Residue+j1,(i1-1)*N_Residue+j1+1]; 
        
        Beta(Ind12(N_loop,:)) = F12;
    end
end

%% F13
Ind13 = nan(N_pair,2);

for i2 = 1:N_Strand
    for j2 = 1:N_Residue-2
        
        N_loop = (i2-1)*N_Residue+j2;
        Ind13(N_loop,:) = [(i2-1)*N_Residue+j2, (i2-1)*N_Residue+j2+2]; 
        
        Beta(Ind13(N_loop,:)) = F13;
    end
end

%% Fab
Ind_ab = nan(Num_Modes - N_Residue,2);

for i3 = 1:length(Ind_ab)
    Ind_ab(i3,:) = [i3,i3+N_Residue];
    Beta(Ind_ab(i3,:)) = Fab;
end

%% Fac
Ind_ac = nan(Num_Modes - 2*N_Residue,2);

for i4 = 1:length(Ind_ac)
    Ind_ac(i4,:) = [i4,i4+2*N_Residue];
    Beta(Ind_ac(i4,:)) = Fac;
end

%% Fpb1
for i5 = 1:N_Strand-1
    for j5 = 2:2:N_Residue
        Indpb1 = j5 + (i5-1)*N_Residue;
        
        Beta(Indpb1,Indpb1 + N_Residue -1) = Fpb1;
    end
end

%% Fpb2
for i6 = 1:N_Strand-1
    for j6 = 1:2:N_Residue-1
        Indpb2 = j6 + (i6-1)*N_Residue;
        
        Beta(Indpb2,Indpb2 + N_Residue +1) = Fpb2;
    end
end

%% Fpb3
for i7 = 1:N_Strand-1
    for j7 = 2:2:N_Residue-1
        Indpb3 = j7 + (i7-1)*N_Residue;
        
        Beta(Indpb3,Indpb3 + N_Residue +1) = Fpb3;
    end
end

%% Fpb4
for i8 = 1:N_Strand-1
    for j8 = 1:2:N_Residue
        Indpb4 = j8 + (i8-1)*N_Residue;
        
        Beta(Indpb4,Indpb4 + N_Residue -1) = Fpb4;
    end
end

%% Symmetrize
Beta = Beta + Beta';





