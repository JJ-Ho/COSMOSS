function R = LabFrameAvg(Symmetry,Dimension)

SaveName = ['R',num2str(Dimension),'_precalculated'];

switch Symmetry
    case 'C1'
        % Initially generate using following equation
        % R = hf_R(0,0,0);
        V_Name = 'C1';
        
    case 'C2'
        % Initially generate using following equation
        %R = (hf_R(0 ,0,0)+...
        %     hf_R(pi,0,0))./2;
        V_Name = 'C2_z';

    case 'C4'
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
        disp(['Does not support ',Symmetry,' Symmetry yet...'])
end

%% Load and rename 
load(SaveName,'-mat',V_Name);
S = whos(V_Name,'-file',SaveName);

if isempty(S)
    disp(['No precalculated <R> ' V_Name ' is avalible in ' SaveName ', abort!'])
    R = sparse(zeros(3^(Dimension)));
    return
else
    R = eval(V_Name);
end

% %% clean up numerical errors
% R(abs(R)<1e-10) = 0;