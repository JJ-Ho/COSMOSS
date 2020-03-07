function  hSurf = PlotAtom(hAx,AtomName,XYZ)

%% Switch atom appearence
FaceAlpha = 1;
EdgeColor = 'none';
switch AtomName
    case 'C'
        FaceColor = [0,0,0];
        
    case 'O'
        FaceColor = [1,0,0];
        
    case 'N'
        FaceColor = [0,0,1];
        
    case 'H'
        FaceColor = [1,1,1];
        
    case 'S'
        FaceColor = [1,1,0];
        
    case 'Label'
        FaceColor = [0,1,1];
        FaceAlpha = 0.3;   
    
    otherwise    
        FaceColor = [1,0,1];
        EdgeColor = [0,1,0];
        disp('There is an undefined atom color in PlotAtom.m...')
end

%% determine readius
if strcmp(AtomName,'Label')
    Scaling = 1/100;
    AtomName = 'C';
else
    Scaling = 1/100/2;
end

% create a fake SD to use the SD_AtomicProperties method
SD = StructureData;
SD.AtomName = {AtomName,AtomName};
SD.XYZ = zeros(2,3);
AP = SD_AtomicProperties(SD);
Radii = AP.Radii(1);
R = Scaling*Radii;

%% draw
[x,y,z] = sphere;
N = size(XYZ,1);

for j=1:N
  hSurf = surf(hAx,R*x+XYZ(j,1),R*y+XYZ(j,2),R*z+XYZ(j,3));
  set(hSurf,...
      'FaceColor',FaceColor,...
      'FaceAlpha',FaceAlpha,...
      'FaceLighting','gouraud',...
      'EdgeColor',EdgeColor) 
end
