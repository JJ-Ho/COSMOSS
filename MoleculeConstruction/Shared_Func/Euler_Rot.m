function Rot = Euler_Rot(XYZ_rotated_frame,xyz_fixed_frame)
% 
%  According to wiki's notation:
%   lowercase x,y,z refers to fixed axes 
%  (original) Uppercase X,Y,Z refers to rotated axes (new)
%  Given any orthonormal coorodinat basis XYZ_Lab_frame, the Eular_Rot 
%  function will gives the rotational matrix to rotate molecular frame
%  (indentity matrix) into Lab frame.
% 
% ------- Version log -----------------------------------------------------
% 
% Ver. 2.1  131028  Modified gamma
% 
% Ver. 2.0  131028  Angles and Rotations matrix are fixed according to
%                   wiki's notation = > active intrinsic rotation.
% 
% Ver. 1.0  130604  Modified from Eular_rot in YCC_Lab folder
% 
% ------------------------------------------------------------------------
% Copyright Jia-Jung Ho, 2013

%% debug
% XYZ_rotated_frame = eye(3);
% xyz_fixed_frame = LAB_Frame;

%% Check if the XYZ_rotated_frame = xyz_fixed_frame
if eq(XYZ_rotated_frame,xyz_fixed_frame)
    Rot = eye(3);
else
    %% Extract Lab frame basis
    v_X = XYZ_rotated_frame(:,1);
    v_Y = XYZ_rotated_frame(:,2);
    v_Z = XYZ_rotated_frame(:,3);

    %% define Molecular frame basis
    % v_x=[1 0 0];
    % v_y=[0 1 0];
    % v_z=[0 0 1];

    v_x = xyz_fixed_frame(:,1)';
    v_y = xyz_fixed_frame(:,2)';
    v_z = xyz_fixed_frame(:,3)';

    % Node defined as intersection of xy and XY planes
    node=cross(v_z,v_Z);
    node=node/norm(cross(v_z,v_Z));

    % Eular angle, a=alpha; b=beta; r=gamma
    %[131028 fixed]
    % alpha=atan2(node*v_y', node*v_x');
    % Y2=cross(node,v_X)/norm(cross(node,v_X));
    % beta=atan2(norm(cross(v_z,v_Z)), dot(v_z,v_Z)); 
    % gamma=atan2(norm(cross(v_X,node)), dot(v_X,node));
    alpha=atan2(dot(node,v_y), dot(node,v_x));
    beta=atan2(dot(v_Z,cross(node,v_z)), dot(v_Z,v_z));
    % gamma=atan2(dot(v_Y,-node), dot(v_X,node));
    gamma=atan2(dot(v_X,cross(v_Z,node)), dot(v_X,node));

    c1=cos(alpha);
    s1=sin(alpha);
    c2=cos(beta);
    s2=sin(beta);
    c3=cos(gamma);
    s3=sin(gamma);

    % The rotation matrix follows Z-X-Z convention of intrinsic rotation  
    % Rot=[ c1*c3-c2*s1*s3  c1*s3+c2*c3*s1  s1*s2
    %      -c3*s1-c1*c2*s3 -s1*s3+c1*c2*c3  c1*s2
    %       s2*s3          -c3*s2           c2    ];

    % 131023 ccording to wiki: Euler_angles
    % Active, intrinsic rotation 
    Rot=[ c1*c3-c2*s1*s3 -c1*s3-c2*c3*s1  s1*s2
          c3*s1+c1*c2*s3  c1*c2*c3-s1*s3 -c1*s2
          s2*s3           c3*s2           c2    ];

end