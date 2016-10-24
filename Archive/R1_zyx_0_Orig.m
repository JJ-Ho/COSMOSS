function M_zyx = R1_zyx_0(D_X,D_Y,D_Z)

c1=cosd(D_X);
s1=sind(D_X);
c2=cosd(D_Y);
s2=-sind(D_Y);
c3=cosd(D_Z);
s3=-sind(D_Z);


% According to wiki: Euler_angles, 
% Tait-Bryan angles, following z-y'-x" 
M_zyx=[ c1*c2, c1*s2*s3 -    c3*s1,    s1*s3 + c1*c3*s2;
        c2*s1,    c1*c3 + s1*s2*s3, c3*s1*s2 -    c1*s3;
          -s2,               c2*s3,               c2*c3];