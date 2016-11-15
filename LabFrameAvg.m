function R = LabFrameAvg(Symmetry,Dimension)

R_name = ['R',num2str(Dimension),'_ZYZ_0'];
hf_R = str2func(R_name);

switch Symmetry
    case 'No'
        R = hf_R(0,0,0);
    case 'C2'
        R = (hf_R(0 ,0,0)+...
             hf_R(pi,0,0))./2;

    case 'C4'
        R = (hf_R(0     ,0,0)+...
             hf_R(pi/2  ,0,0)+...
             hf_R(pi    ,0,0)+...
             hf_R(pi/2*3,0,0))./4;
    case 'Isotopic'
end

% clean up numerical errors
R(abs(R)<1e-10) = 0;