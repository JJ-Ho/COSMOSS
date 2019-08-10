function FormattedPath = TableFromatPathways(PathType,PathInd)
%% Debug
% clear all
% PathType = [1,2,3,4,6];
% PathInd = reshape(1:20,5,4);

%% Pre allocate cell
N_Path = size(PathInd,1);
I1 = cell(N_Path,1);
I2 = cell(N_Path,1);
I3 = cell(N_Path,1);
I4 = cell(N_Path,1);

%% Detemine PathType logical index
IR1  = eq(PathType,1);
IR2  = eq(PathType,2);
IR3  = eq(PathType,3);
INR1 = eq(PathType,4);
INR2 = eq(PathType,5);
INR3 = eq(PathType,6);

%% Number format
FS = '%2u';
FZ = ' 0';

%%  R1, GB
%  |b0|
%  |00|
%  |0a|
a = PathInd(IR1,4);
b = PathInd(IR1,2);

I1(IR1) = strcat(FZ,{','},num2str(a,FS));
I3(IR1) = strcat(num2str(b,FS),{','},FZ);

[I2{IR1}] = deal([FZ,',',FZ]);
[I4{IR1}] = deal([FZ,',',FZ]);

%% R2, SE
% |b0| 
% |ba|
% |0a|
a = PathInd(IR2,4);
b = PathInd(IR2,3);

I1(IR2) = strcat(FZ,{','},num2str(a,FS));
I2(IR2) = strcat(num2str(b,FS),{','},num2str(a,FS));
I3(IR2) = strcat(num2str(b,FS),{','},FZ);

[I4{IR2}] = deal([FZ,',',FZ]);        
        
%% R3, EA
% |xa| 
% |ba|
% |0a|
a = PathInd(IR3,4);
b = PathInd(IR3,3);
x = PathInd(IR3,2);

I1(IR3) = strcat(FZ,{','},num2str(a,FS));
I2(IR3) = strcat(num2str(b,FS),{','},num2str(a,FS));

I3(IR3) = strcat('<html><pre style="color:red ;">', ...
           strcat(num2str(x,FS),{','},num2str(a,FS)), ...
           '</pre></html>');

I4(IR3) = strcat(num2str(a,FS),{','},num2str(a,FS));

%% NR1, GB
%  |b0|
%  |00|
%  |a0|
a = PathInd(INR1,4);
b = PathInd(INR1,2);

I1(INR1) = strcat(num2str(a,FS),{','},FZ);
I3(INR1) = strcat(num2str(b,FS),{','},FZ);

[I2{INR1}] = deal([FZ,',',FZ]);
[I4{INR1}] = deal([FZ,',',FZ]);

%% NR2, SE
% |a0| 
% |ab|
% |a0|
a = PathInd(INR2,4);
b = PathInd(INR2,3);

I1(INR2) = strcat(num2str(a,FS),{','},FZ);
I2(INR2) = strcat(num2str(a,FS),{','},num2str(b,FS));
I3(INR2) = strcat(num2str(a,FS),{','},FZ);

[I4{INR2}] = deal([FZ,',',FZ]);     

%% NR3, EA
% |xb| 
% |ab|
% |a0|
a = PathInd(INR3,4);
b = PathInd(INR3,3);
x = PathInd(INR3,2);

I1(INR3) = strcat(num2str(a,FS),{','},FZ);
I2(INR3) = strcat(num2str(a,FS),{','},num2str(b,FS));

I3(INR3) = strcat('<html><pre style="color:red ;">', ...
           strcat(num2str(x,FS),{','},num2str(b,FS)), ...
           '</pre></html>');

I4(INR3) = strcat(num2str(b,FS),{','},num2str(b,FS));


FormattedPath = [I4,I3,I2,I1];