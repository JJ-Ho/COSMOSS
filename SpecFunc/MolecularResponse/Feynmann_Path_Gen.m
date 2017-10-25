function [SpectraGrid,Beta,IntensityGrid] = Feynmann_Path_Gen(SpecType,Pathways,Data_2D,SparseMax,MEM_CutOff)
%% debug
% SpecType = '2DSFG';
% Pathways = 'NR3';
% 
% S = Data_COSMOSS.Structure;
% I = Data_COSMOSS.hCOSMOSS.Parse_GUI;
% H = ExcitonH(S,I,'TwoEx');    
% 
% Mu    = MuAlphaGen(S,H,'Mode','Mu');
% Alpha = MuAlphaGen(S,H,'Mode','Alpha');
% 
% Data_2D.H = H;
% Data_2D.Mu = Mu;
% Data_2D.Alpha = Alpha;
% Data_2D.PCutOff = 1e-5;
% SparseMax = 2000;
% MEM_CutOff = 3e-1; %[GB]

%% Define common variables
H       = Data_2D.H;
Mu      = Data_2D.Mu;
PCutOff = Data_2D.PCutOff;
EJLR    = Data_2D.EJLR;

F1    = H.Sort_Ex_F1;
F2    = H.Sort_Ex_F2;
M01   = Mu.M_Ex_01;
M01_N = Mu.M_Ex_01_N;

%% Decide Frequency bin size
F1 = round(F1); % 1cm-1
F2 = round(F2); 

%% Assign variables according to Pathway and Spectral type
N1 = length(F1);
N2 = length(F2);

% Prepare interaction indexes for the correspinding pathway
switch Pathways
    case {'R1','R2','NR1','NR2'}
        [Ib,Ia] = ndgrid(1:N1,1:N1);
        Ia = Ia(:); % 1ex idx
        Ib = Ib(:); % 1ex idx
        
    case {'R3','NR3'}
        [Kx,Kb,Ka] = ndgrid(1:N2,1:N1,1:N1);
        Ka = Ka(:); % 1ex idx
        Kb = Kb(:); % 1ex idx
        Kx = Kx(:); % 2ex idx
        
        % linearlized  1ex -> 2ex transition matrix Index
        % Note: The length of Iax/Ibx = N1^3*(N1+1)/2
        Iax = sub2ind([N1,N2],Ka,Kx);
        Ibx = sub2ind([N1,N2],Kb,Kx);
        
        % Define the linearlized  1ex -> 2ex transition dipole matrix 
        M12 = reshape(Mu.M_Ex_12,[],3);
        M12_N = Mu.M_Ex_12_N;
        
    otherwise
        disp(['Pathway type: ',Pathways,' is not correct...'])
        return
end

% Assign the transition index according to the selected pathway
switch Pathways
    case {'R1','NR1'}
        I1 = Ia;
        I2 = Ia;
        I3 = Ib;
        I4 = Ib;
    case 'R2'
        I1 = Ia;
        I2 = Ib;
        I3 = Ia;
        I4 = Ib;
    case 'NR2'
        I1 = Ia;
        I2 = Ib;
        I3 = Ib;
        I4 = Ia;        
    case 'R3'
        I1 = Ka;
        I2 = Kb;
        I3 = Ibx;
        I4 = Iax; 
    case 'NR3'
        I1 = Ka;
        I2 = Kb;
        I3 = Iax;
        I4 = Ibx; 
end

% Give tensor indexes and define the last interaction with respect to the
% corresponding Spectral type
switch SpecType
    case '2DSFG'
        [Jd,Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3,1:9);

        L01   = Data_2D.Alpha.M_Ex_01;
        L01_N = Data_2D.Alpha.M_Ex_01_N;
        L12   = reshape(Data_2D.Alpha.M_Ex_12,[],9); % linearlized  1ex -> 2ex polarizability matrix
        L12_N = Data_2D.Alpha.M_Ex_12_N;
        
    case '2DIR'
        [Jd,Jc,Jb,Ja] = ndgrid(1:3,1:3,1:3,1:3);
        
        L01   = M01;
        L01_N = M01_N;
        L12   = reshape(Mu.M_Ex_12,[],3); % linearlized  1ex -> 2ex transition dipole matrix 
        L12_N = Mu.M_Ex_12_N;
        
    otherwise
        disp(['Spectral type: ',SpecType,' is not correct...'])
        return
end
Ja = Ja(:);
Jb = Jb(:);
Jc = Jc(:);
Jd = Jd(:);

L_Response = length(Ja);

% Decide Transition vectors for the corresponding pathway type
switch Pathways
    case {'R1','R2','NR1','NR2'}
        T1 = M01; T1N = M01_N;
        T2 = M01; T2N = M01_N;
        T3 = M01; T3N = M01_N;
        T4 = L01; T4N = L01_N;
        
    case {'R3','NR3'}       
        T1 = M01; T1N = M01_N;
        T2 = M01; T2N = M01_N;
        T3 = M12; T3N = M12_N;
        T4 = L12; T4N = L12_N;
end

%% Reduce mode numbers base on the Response intensity
% Norm_T1 = sqrt(sum(T1.^2,2));
% Norm_T2 = sqrt(sum(T2.^2,2));
% Norm_T3 = sqrt(sum(T3.^2,2));
% Norm_T4 = sqrt(sum(T4.^2,2));
% Intentisy = Norm_T4(I4).*Norm_T3(I3).*Norm_T2(I2).*Norm_T1(I1);

Intentisy = T4N(I4).*T3N(I3).*T2N(I2).*T1N(I1);
Int_Max = max(Intentisy);

% Apply CutOff
Reduced_Ind = Intentisy  >= (PCutOff *  Int_Max);
I1 = I1(Reduced_Ind);
I2 = I2(Reduced_Ind);
I3 = I3(Reduced_Ind);
I4 = I4(Reduced_Ind);

%% Generate [Pump,Probe] Frequencies Coordinates
switch Pathways
    case {'R1','R2','NR1','NR2'}
        F_sub = [ F1(I1), F1(I4)]; 
    case 'R3'
        Kx = Kx(Reduced_Ind); % note: non-linearlized "1ex -> 2ex" index for F2
        F_sub = [ F1(I1), F2(Kx) - F1(I1)];
    case 'NR3'
        Kx = Kx(Reduced_Ind); % note: non-linearlized "1ex -> 2ex" index for F2
        F_sub = [ F1(I1), F2(Kx) - F1(I2)];
end

% Linearlize Frequency from 2D to 1D Grid
% This way, I can generate sparsed Beta matrix
F_ind = sub2ind([SparseMax,SparseMax],F_sub(:,1),F_sub(:,2)); 

%% Memory cutoff, estimate the largest array (Beta)'s size and break it down to several for-loop
Ele_Max = round(MEM_CutOff/(L_Response * 8 / 1e9)) + 1; % max number of elements to reach MEM_CufOff

%% Calculate molecular response tensor
Beta = sparse(SparseMax^2,L_Response); % Preallocate sparse matrix container
N_Path = numel(I1);

if N_Path > Ele_Max
    % Add NaNs so I can reshape indexes
    Padding_L = Ele_Max - mod(N_Path,Ele_Max);
    Padding_NaN  = nan(Padding_L,1);
    
    F1D_Loop = reshape([F_ind; Padding_NaN],Ele_Max,[]);
    I1 = reshape([I1; Padding_NaN],Ele_Max,[]);
    I2 = reshape([I2; Padding_NaN],Ele_Max,[]);
    I3 = reshape([I3; Padding_NaN],Ele_Max,[]);
    I4 = reshape([I4; Padding_NaN],Ele_Max,[]);

    Loop_N = size(I1,2);
    % display how many loop is going to be done
    MEM = N_Path * L_Response * 8 / 1e9; % roughly unit in GB
    disp('--------------------------------------')
    disp(['Memory cut-off set at ', num2str(MEM_CutOff), 'GB...'])
    disp(['Memory size of ',Pathways,' is about ', sprintf('%4.1f',MEM), 'GB...'])
    disp(['Will run ',Pathways,' with ', num2str(Loop_N), ' Loops to reduce memory load...'])
    disp('--------------------------------------')
    
    Response = zeros(Ele_Max,Loop_N-1);
    for LN = 1:Loop_N -1         
        Beta_Raw = T4(I4(:,LN),Ja).*T3(I3(:,LN),Jb).*T2(I2(:,LN),Jc).*T1(I1(:,LN),Jd);
        Beta(F1D_Loop(:,LN),:) = Beta_Raw;
        Response(:,LN)         = Beta_Raw*EJLR';
    end    
    Response = Response(:);
    
    % deall with the last loop
    Last_Ind = (1:mod(N_Path,Ele_Max))';
    F1D_End  = F1D_Loop(Last_Ind,Loop_N);
    I1_End   = I1(Last_Ind,Loop_N);
    I2_End   = I2(Last_Ind,Loop_N);
    I3_End   = I3(Last_Ind,Loop_N);
    I4_End   = I4(Last_Ind,Loop_N);
    
    Beta_Raw_End = T4(I4_End,Ja).*T3(I3_End,Jb).*T2(I2_End,Jc).*T1(I1_End,Jd);
    Beta(F1D_End,:) = Beta_Raw_End;
    Response_End    = Beta_Raw_End*EJLR';
    
    Response = [ Response; Response_End];    
else
    Beta_Raw = T4(I4,Ja).*T3(I3,Jb).*T2(I2,Jc).*T1(I1,Jd);
    Beta(F_ind,:) = Beta_Raw;
    Response      = Beta_Raw*EJLR';
end

%% Deal with Other outputs
SpectraGrid   = sparse(F_sub(:,1),F_sub(:,2),Response,SparseMax,SparseMax);

% IntensityGrid = reshape((sum((Beta).^2,2).^(1/2)),SparseMax,SparseMax);

% regenerate intensity after cutoff
IntCutOff = T4N(I4).*T3N(I3).*T2N(I2).*T1N(I1);
IntensityGrid = sparse(F_sub(:,1),F_sub(:,2),IntCutOff,SparseMax,SparseMax);

