function Output = ConstructGrid(GUI_Inputs)
%% ConstructGrid

% ------- Version log -----------------------------------------------------
% 
% Ver. 1.1  150426  Copy from 2DSFG_Ester project
%                   Update debug part
%                   Implement inputparser
% 
% Ver. 1.0  140903  Modify from GetAcid.m 
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2014

%% Debug input part
% clear all
% GUI_Inputs.debug = 1;

%% Inputs parser
GUI_Inputs_C      = fieldnames(GUI_Inputs);
GUI_Inputs_C(:,2) = struct2cell(GUI_Inputs);
GUI_Inputs_C      = GUI_Inputs_C';

INPUT = inputParser;
INPUT.KeepUnmatched = 1;

% Default values
defaultAng_Phi     = 0;
defaultAng_Psi     = 0;
defaultAng_Theta   = 0;
defaultDelta_Phi   = 0;
defaultDelta_Psi   = 0;
defaultDelta_Theta = 0;
defaultVec_1       = [7,0,0];
defaultVec_2       = [0,4,0];
defaultN_1         = 2;
defaultN_2         = 3;
defaultNLFreq      = 1720;
defaultMute_Ind    = 'Non';

% add Optional inputs / Parameters
addOptional(INPUT,'Ang_Phi'    ,defaultAng_Phi    );
addOptional(INPUT,'Ang_Psi'    ,defaultAng_Psi    );
addOptional(INPUT,'Ang_Theta'  ,defaultAng_Theta  );
addOptional(INPUT,'Delta_Phi'  ,defaultDelta_Phi  );
addOptional(INPUT,'Delta_Psi'  ,defaultDelta_Psi  );
addOptional(INPUT,'Delta_Theta',defaultDelta_Theta);
addOptional(INPUT,'Vec_1'      ,defaultVec_1      );
addOptional(INPUT,'Vec_2'      ,defaultVec_2      );
addOptional(INPUT,'N_1'        ,defaultN_1        );
addOptional(INPUT,'N_2'        ,defaultN_2        );
addOptional(INPUT,'NLFreq'     ,defaultNLFreq     );
addOptional(INPUT,'Mute_Ind'   ,defaultMute_Ind   );

parse(INPUT,GUI_Inputs_C{:});

% Re-assign variable names
Ang_Phi     = INPUT.Results.Ang_Phi    ;
Ang_Psi     = INPUT.Results.Ang_Psi    ;
Ang_Theta   = INPUT.Results.Ang_Theta  ;
Delta_Phi   = INPUT.Results.Delta_Phi  ;
Delta_Psi   = INPUT.Results.Delta_Psi  ;
Delta_Theta = INPUT.Results.Delta_Theta;
Vec_1       = INPUT.Results.Vec_1      ;
Vec_2       = INPUT.Results.Vec_2      ;
N_1         = INPUT.Results.N_1        ;
N_2         = INPUT.Results.N_2        ;
NLFreq      = INPUT.Results.NLFreq     ;
Mute_Ind = INPUT.Results.Mute_Ind;

%% Read G09 structure and reassign variables
    
RR = [Ang_Phi,Ang_Psi,Ang_Theta];
RR = RR./180*pi; % turn to radius unit



P = ReadG09Input('140903_MMB.txt',RR);
XYZ       = P.XYZ;
Atom_Num  = P.Atom_Num;
mu_Mol    = P.TDV;
alpha_Mol = P.Raman;

% % shift monomer to origin
% XYZ = bsxfun(@minus,XYZ,sum(XYZ,1)/Atom_Num);
%% Define usefule constants

% define 2D grid vector
% Vec_1 = a.*[1,0,0];
% Vec_2 = b.*[0,1,0];

Num_Modes = N_1*N_2;

Phi_Fluc   = Delta_Phi*randn(Num_Modes,1)  ./180*pi;
Psi_Fluc   = Delta_Psi*randn(Num_Modes,1)  ./180*pi;
Theta_Fluc = Delta_Theta*randn(Num_Modes,1)./180*pi;

%% Creating Translation Copy
XYZ_Grid_M = zeros(Atom_Num,3,N_1,N_2);

mu_Sim    = zeros(Num_Modes,3);
alpha_Sim = zeros(Num_Modes,3,3);


for j = 1:N_2
    for i = 1:N_1

        if any(Mute_Ind == i+(j-1)*N_1)
            % Mute some molecule in the grid
            XYZ_Grid_M(:,:,i,j)      = zeros(Atom_Num,3);
            mu_Sim(i+(j-1)*N_1,:,:)    = zeros(3,1);
            alpha_Sim(i+(j-1)*N_1,:,:) = zeros(3,3);
        else     
            % Add random orientation fluctuation
            R_Fluc = R1_ZYZ_0(Phi_Fluc(i+(j-1)*N_1),Psi_Fluc(i+(j-1)*N_1),Theta_Fluc(i+(j-1)*N_1));

            % XYZ
            XYZ_R  = (R_Fluc*XYZ')'; % apply random fluctuation 
            TransV = (i-1)*Vec_1 + (j-1)*Vec_2;
            XYZ_Grid_M(:,:,i,j) = bsxfun(@plus,XYZ_R,TransV);

            % Mu & Alpha
               mu_Sim(i+(j-1)*N_1,:,:) = (R_Fluc*mu_Mol')';
            alpha_Sim(i+(j-1)*N_1,:,:) =  R_Fluc*alpha_Mol*R_Fluc';
        end
    end
end

XYZ_Grid = reshape(permute(XYZ_Grid_M,[1,3,4,2]),[],3);


%% Make debug figure

% Atom_Num_NM = Atom_Num*N_1*N_2;
% 
% [n,m] = ndgrid(1:Atom_Num_NM);
% T1=XYZ_Grid(m(:),:)-XYZ_Grid(n(:),:);
% T2=sqrt(sum(T1.^2,2));
% Distance_matrix=reshape(T2,Atom_Num_NM,Atom_Num_NM);
% lower=tril(Distance_matrix,-1);
% [a,b]=find(lower<1.6 & lower>0);
% C_index=[a,b];
% Conn_Grid=false(Atom_Num_NM);
% Conn_Grid(C_index(:,1)+(C_index(:,2)-1)*Atom_Num_NM)=true;
% 
% for k = 1:Num_Modes
%     Conn_Grid(13+(k-1)*Atom_Num,3+(k-1)*Atom_Num) = 1; % Add S-C connection for Ester
% end
% 
% Conn_Grid=Conn_Grid|Conn_Grid';
% figure
% gplot3(Conn_Grid,XYZ_Grid)
% axis equal

%% Create Translational copy of Center
C_Atom_Ind = 7;

Center_M = XYZ_Grid_M(C_Atom_Ind,:,:,:);
Center = reshape(permute(Center_M,[1,3,4,2]),[],3);

%% Calculate the transition dipoles (mu) and Raman tensors (alpha) in the Lab frame
% 
% mu_Mol    = P.TDV;
% alpha_Mol = P.Raman;
% 
% mu_Sim    = zeros(Num_Modes,3);
% alpha_Sim = zeros(Num_Modes,3,3);
% 
% for ii=1:Num_Modes
%     mu_Sim(ii,:,:)    = (mu_Mol);
%     alpha_Sim(ii,:,:) = alpha_Mol;
% end 

%% Define Mode frequency and anharmonicity

% if isstruct(handles)
%     NLFreq = str2double(get(handles.E_NLFreq  ,'String'));
%     Anharm = str2double(get(handles.E_Anharm  ,'String'));
% else
%     NLFreq = P.Freq.Fundamental;
%     Anharm = 2*P.Freq.Fundamental - P.Freq.Overtone;
% end
%     
% 
Freq   = ones(Num_Modes,1)* NLFreq;
Anharm = ones(Num_Modes,1)*(2*P.Freq.Fundamental - P.Freq.Overtone);

%% AtomSerNo

AtomSerNo = zeros(Num_Modes,3);
for jj = 1:Num_Modes
    Shift_Ind = (jj-1)*Atom_Num;
    AtomSerNo(jj,:) = [7+Shift_Ind,10+Shift_Ind,8+Shift_Ind];
end

%% Deal with files name 
Grid_FilesName = [ 'Ester_Grid_V1' num2str(N_1) 'V2' num2str(N_2)];

%% Output Structure
Output.center       = Center;
Output.freq         = Freq;
Output.anharm       = Anharm;
Output.mu           = mu_Sim;
Output.alpha        = reshape(alpha_Sim,[Num_Modes,9]); % non-reduced alpha "vector"
Output.alpha_matrix = alpha_Sim;
Output.AtomSerNo    = AtomSerNo;
Output.Num_Modes    = Num_Modes;
Output.XYZ          = XYZ_Grid;
Output.FilesName    = Grid_FilesName;
Output.Monomer      = P;