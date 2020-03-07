function Conn = SD_Connectivity(obj_SD)
N   = obj_SD.NAtoms;
XYZ = obj_SD.XYZ;
AP  = SD_AtomicProperties(obj_SD);
RC  = AP.Radii;

%% generate connectivity matrix ndgrid version
% Define bonded if dist. < (radii A + radii B)* BL_Scale
BL_Scale = 1.5;
% generate index for substracting
[n,m] = ndgrid(1:N);
% Generate bond length matrix, /100 frm pm to Ang, * scaling facotr
B = reshape((RC(m(:)) + RC(n(:)) )./100.*BL_Scale,N,N);
% substract atom's position according to index generate above
T1 = XYZ(m(:),:)-XYZ(n(:),:);
% Get atom distance
T2 = sqrt(sum(T1.^2,2));
% reshape a (anum^2)x1 matrix to anum x anum matrix
Distance_matrix = reshape(T2,N,N);
% Check Bondlength
Bonded = Distance_matrix < B;
% grep lower triangular part in order to avoid over counting
lower = tril(Bonded,-1);
% Transform into bonding index metrix
[a,b] = find(lower>0);
C_index = [a,b];
Conn = false(N);
Conn(C_index(:,1)+(C_index(:,2)-1)*N) = true;
Conn = Conn|Conn';

