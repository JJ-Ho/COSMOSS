clear m;
close all

t=1;
t_max = 500;
SampleRate = 20; % Hz
XArray = 1:t_max;
O = zeros(length(XArray),3);

hF = figure;
L = 40;
hAx = axes(hF,'XLim',[-L,L],'YLim',[-L,L],'ZLim',[-L,L]);
view([131,27])
rotate3d on

m = mobiledev;


m.SampleRate = SampleRate;
m.OrientationSensorEnabled = 1;
m.Logging = 1;

pause(1)

while(t<t_max)
    pause(1/SampleRate)
    
    [Read,T_Read] = orientlog(m);
    O(t,:) = Read(end,:);
    
    %% Prep molecule
    XYZ_0 = Structure.XYZ;
    XYZ_0 = bsxfun(@minus,XYZ_0,sum(XYZ_0,1)./size(XYZ_0,1));
%     XYZ_0 = [1,0,0;
%              0,1,0;
%              0,0,1 ];

    % Orientation = Orientation/180*pi; % turn to radius unit
    DegreeX = O(t,1);
    DegreeY = O(t,2);
    DegreeZ = O(t,3);

    %M_X = Rx(RadiusX);
    %M_Y = Rx(RadiusY);
    %M_Z = Rx(RadiusZ);
    %M = R1_ZYZ_0(RadiusX,RadiusY,RadiusZ);
    
    %% z-y'-x"
    c1=cosd(DegreeX);
    s1=sind(DegreeX);
    c2=cosd(DegreeY);
    s2=-sind(DegreeY);
    c3=cosd(DegreeZ);
    s3=-sind(DegreeZ);


    % According to wiki: Euler_angles, 
    % Tait-Bryan angles, following z-y'-x" 
    M_zyx=[ c1*c2, c1*s2*s3 -    c3*s1,    s1*s3 + c1*c3*s2;
            c2*s1,    c1*c3 + s1*s2*s3, c3*s1*s2 -    c1*s3;
              -s2,               c2*s3,               c2*c3];
   
    XYZ_1 = (M_zyx*XYZ_0')';
    
    Carbon_Pos   = XYZ_1(Structure.AtomSerNo(:,1),:);
    Oxygen_Pos   = XYZ_1(Structure.AtomSerNo(:,2),:);
    Nitrogen_Pos = XYZ_1(Structure.AtomSerNo(:,3),:);
    
    %% make figure
    cla
    hold on
    Conn = Connectivity(XYZ_1);
    gplot3(Conn,XYZ_1);
    
    plot3(Carbon_Pos(:,1)  ,Carbon_Pos(:,2)  ,Carbon_Pos(:,3)  ,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerSize',5)
    plot3(Oxygen_Pos(:,1)  ,Oxygen_Pos(:,2)  ,Oxygen_Pos(:,3)  ,'LineStyle','none','Marker','o','MarkerFaceColor','r','MarkerSize',5)
    plot3(Nitrogen_Pos(:,1),Nitrogen_Pos(:,2),Nitrogen_Pos(:,3),'LineStyle','none','Marker','o','MarkerFaceColor','b','MarkerSize',5)
    
    hold off
    %plot(hAx,XArray(1:t),O(1:t,1),XArray(1:t),O(1:t,2),XArray(1:t),O(1:t,3))
    
    %PlotRotMolFrame(hAx,XYZ_0,XYZ_1,[0,0,0])
    
    grid on
    hAx.Box = 'on';
    axis equal
    hAx.XLim = [-L,L];
    hAx.YLim = [-L,L];
    hAx.ZLim = [-L,L];
    drawnow
    
    t = t+1
end

m.Logging = 0;
m.OrientationSensorEnabled = 0;
clear m;