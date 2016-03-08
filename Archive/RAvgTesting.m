% 131014 All the rotational notation are fixed. Ref wiki page


tic
syms Phi Psi Theta

% Follow wiki's notation http://en.wikipedia.org/wiki/Rotation_matrix

RZ_Psi   = [ cos(Psi), sin(Psi), 0; 
            -sin(Psi), cos(Psi), 0; 
                   0,         0, 1];
        
% RX_Theta = [1,          0,          0; 
%             0, cos(Theta), sin(Theta); 
%             0,-sin(Theta), cos(Theta)]; 

RY_Theta = [cos(Theta), 0, sin(Theta);
            0,          1,          0;   
            0,-sin(Theta), cos(Theta)]; 

RZ_Phi   = [ cos(Phi), sin(Phi), 0;
            -sin(Phi), cos(Phi), 0; 
                   0,         0, 1];

% => according to active intrinsic transformation, rotate molecule
R = RZ_Psi*RY_Theta*RZ_Phi; 
disp('R matrix generated...')

toc
%%
% RR = outer(R,R,0);
% RR = permute(RR,[1,3,2,4]);
% disp('RR matrix generated...')
% RRR = outer(RR,R,0);
% RRR = permute(RRR,[1,2,5,4,3,6]);
% RRR = permute(RRR,[1,2,3,5,4,6]);
% disp('RRR matrix generated...')
% 
% RRRR = outer(RRR,R,0);
% RRRR = permute(RRRR,[1,2,3,7,5,6,4,8]);
% RRRR = permute(RRRR,[1,2,3,4,5,7,6,8]);
% RRRR = permute(RRRR,[1,2,3,4,6,5,7,8]);
% disp('RRRR matrix generated...')
% 
% RRRRR = outer(RRRR,R,0);
% RRRRR = permute(RRRRR,[1,2,3,4,9,6,7,8,5,10]);
% RRRRR = permute(RRRRR,[1,2,3,4,5,6,7,9,8,10]);
% RRRRR = permute(RRRRR,[1,2,3,4,5,6,8,7,9,10]);
% RRRRR = permute(RRRRR,[1,2,3,4,5,7,6,8,9,10]);
% disp('RRRRR matrix generated...')
% 
% RR = reshape(RR,[],3^2);
% RRR = reshape(RRR,[],3^3);
% RRRR = reshape(RRRR,[],3^4);
% RRRRR = reshape(RRRRR,[],3^5);
% disp('All matrix reshaped ...')
% 
% disp('All symbolic matrixies generated...')




%% R1

tic
R1_ZYZ_1 = sym('E',[1,3^2]); 

for i=1:3^2
    tic
    R1_ZYZ_1(i) = int(R(i),Phi,0,2*pi)./(2*pi);
    
    str = sprintf('R1_ZYZ_1 %d-th run is finished...', i);
    disp(str);
    toc
end
R1_ZYZ_1 = reshape(R1_ZYZ_1,3,3);

% % Save all variables
% C=clock;
% Year = num2str(C(1)-2000); % shift 20XX to XX for shorter filename
% 
% if C(2)<10
%     Month = ['0' num2str(C(2))];
% else
%     Month = num2str(C(2));
% end
% 
% if C(3)<10
%     Date = ['0' num2str(C(3))];
% else
%     Date = num2str(C(3));
% end
% 
% Time = num2str(C(4:6));
% str = [Year  Month Date '_XZX_R5_1' ];
% save(str)
% 
% disp('Integration part finished, please check any Warning in the Output file')


%% R5 saved into function
% tic
% matlabFunction(R1_ZYZ_1,'file','R5_ZXZ_1');
% disp('R5_1 saved as Matlab function!')
% toc

% matlabpool close