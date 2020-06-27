function [R,M,R_List,M_List] = LabFrameAvg(Avg_Rot,Avg_Mirror,N_Interactions)
%% List of options
R_List = {'Isotropic','C1','C2_z','C4_z (<Phi>)'};
M_List = {'No','Sigma v'};

R = 'Not Assigned';
M = 'Not Assigned';
%% Load precalculated Rotational average and rename 
switch Avg_Rot
    case 'List'
        V_Name = 'None';
        
    case 'C1'
        % Initially generate using following equation
        % R = hf_R(0,0,0);
        V_Name = 'C1';
        
    case 'C2_z'
        % Initially generate using following equation
        %R = (hf_R(0 ,0,0)+...
        %     hf_R(pi,0,0))./2;
        V_Name = 'C2_z';

    case 'C4_z (<Phi>)'
        % Initially generate using following equation
        % R = (hf_R(0     ,0,0)+...
        %      hf_R(pi/2  ,0,0)+...
        %      hf_R(pi    ,0,0)+...
        %      hf_R(pi/2*3,0,0))./4;
        % To save time, I load the precalculated matrix instead
        V_Name = 'C4_z';
        
    case 'Isotropic'
        % Initially generate by integration of 
        % Rn_ZYZ_0(phi,psi,theta)*sin(theta)d(phi)d(psi)d(theta)
        V_Name = 'C_iso';
        
    otherwise
        disp(['Does not support ',Avg_Rot,' Symmetry yet...'])
end

SaveName = ['R',num2str(N_Interactions),'_precalculated'];
if ~strcmp(V_Name,'None')
    load(SaveName,'-mat',V_Name);
    S = whos(V_Name,'-file',SaveName);

    if isempty(S)
        disp(['No precalculated <R> ' V_Name ' is avalible in ' SaveName ', abort!'])
        R = sparse(zeros(3^(N_Interactions)));
        return
    else
        R = eval(V_Name);
    end
end

%% Decide Mirror planes
switch Avg_Mirror
    case 'No' % no mirror plane
        V = [1;1;1];
        switch N_Interactions
            case 2
                M = kron(V,V);
            case 3
                M = kron(kron(V,V),V);
            case 5
                M = kron(kron(kron(kron(V,V),V),V),V);
            otherwise
                M = ones(3^(N_Interactions));
        end

    case 'Sigma v' % sigma v, X=-X, Y=-Y
        V1 = [-1; 1;1];
        V2 = [ 1;-1;1];
        switch N_Interactions
            case 2
                Sigma_X = kron(V1,V1);
                Sigma_Y = kron(V2,V2);
            case 3
                Sigma_X = kron(kron(V1,V1),V1);
                Sigma_Y = kron(kron(V2,V2),V2);
            case 5
                Sigma_X = kron(kron(kron(kron(V1,V1),V1),V1),V1);
                Sigma_Y = kron(kron(kron(kron(V2,V2),V2),V2),V2);
            otherwise
                Sigma_X = ones(3^(N_Interactions));
                Sigma_Y = ones(3^(N_Interactions));
        end
        Sigma_X(eq(Sigma_X,-1)) = 0;
        Sigma_Y(eq(Sigma_Y,-1)) = 0;
        
        M = and(Sigma_X,Sigma_Y);
end