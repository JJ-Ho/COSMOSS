function R = Rx(A)
% note positive angle means clockwise rotation

    R = [1,      0,      0; 
         0, cos(A), sin(A); 
         0,-sin(A), cos(A)]; 
end