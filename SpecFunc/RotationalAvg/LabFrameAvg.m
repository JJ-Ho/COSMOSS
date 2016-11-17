function R = LabFrameAvg(Symmetry,Dimension)

R_name = ['R',num2str(Dimension),'_ZYZ_0'];
hf_R = str2func(R_name);

switch Symmetry
    case 'C1'
        R = hf_R(0,0,0);
    case 'C2'
        R = (hf_R(0 ,0,0)+...
             hf_R(pi,0,0))./2;

    case 'C4'
        % Initially generate using following equation
        % R = (hf_R(0     ,0,0)+...
        %      hf_R(pi/2  ,0,0)+...
        %      hf_R(pi    ,0,0)+...
        %      hf_R(pi/2*3,0,0))./4;
        % To save time, I load the pre calculated matrix instead
        SaveName = ['R',num2str(Dimension),'_C4_z'];
        load(SaveName);
        R = C4_z;
        
    case 'Isotropic'
        % Initially generate using following equation
        % C4_z = (hf_R(0     ,0,0)+...
        %         hf_R(pi/2  ,0,0)+...
        %         hf_R(pi    ,0,0)+...
        %         hf_R(pi/2*3,0,0))./4;
        %     
        % C4_x = (hf_R(0,0     ,0)+...
        %         hf_R(0,pi/2  ,0)+...
        %         hf_R(0,pi    ,0)+...
        %         hf_R(0,pi/2*3,0))./4;   
        % 
        % R = sparse(C4_x) * sparse(C4_z);
        % To save time, I load the pre calculated matrix instead
        SaveName = ['R',num2str(Dimension),'_Iso'];
        load(SaveName);
        R = Iso;
        
    otherwise
        disp(['Does not support ',Symmetry,'Symmetry...'])
end

% clean up numerical errors
R(abs(R)<1e-10) = 0;