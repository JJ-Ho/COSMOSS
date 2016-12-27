%% Geenration of E field polarization tensor
% Sig ~ EJLR(Beta)
% 
% This is the E part of the total signal. Given polarizer angle of each 
% incident beams: O_Vis,O_Probe,O_Pump2,O_Pump1, one can get a matrix-lized
% tensor that can apply on Chi_J (Chi_PPPPP, Chi_SPPPP,... ,Chi_SSSSS) and 
% get the Jones vector of the 2DSFG Signal.
% 
% This generation script is written to generate the matrix-lized tensor for
% both 1DSFG and 2DSFG.
% 
% 130911 JJHo
SaveAsFunc = 'y';
syms O_Sig O_Vi O_Pu1 O_Pu2 O_Pr

E_Sig = [cos(O_Sig), sin(O_Sig)];
E_Vis = [cos(O_Vi) , sin(O_Vi) ];
E_Pr  = [cos(O_Pr) , sin(O_Pr) ];
E_Pu2 = [cos(O_Pu2), sin(O_Pu2)];
E_Pu1 = [cos(O_Pu1), sin(O_Pu1)];

E1 =                     E_Sig;
E2 =                kron(E_Sig,E_Vis);
E3 =           kron(kron(E_Sig,E_Vis),E_Pu1);
E4 =      kron(kron(kron(E_Sig,E_Pr ),E_Pu2),E_Pu1);
E5 = kron(kron(kron(kron(E_Sig,E_Vis),E_Pr),E_Pu2),E_Pu1);

%% 
if strcmp(SaveAsFunc,'y')   
    matlabFunction(E1,'file','EPolar1');
    matlabFunction(E2,'file','EPolar2');
    matlabFunction(E3,'file','EPolar3');
    matlabFunction(E4,'file','EPolar4');
    matlabFunction(E5,'file','EPolar5');
end
