function SaveFigures(hF,PathName,FileName)

% check if the folder exixt
if ~exist(PathName,'dir')
    mkdir(PathName)
end

SaveName = [PathName,'/',FileName];
savefig(hF,[SaveName,'.fig'])
saveas (hF,SaveName,'png')

hAx = findobj(hF,'Type','Axes');

NAx = length(hAx);
for i = 1:NAx
    hAx(i).XGrid           = 'off';
    hAx(i).YGrid           = 'off';
    hAx(i).XMinorGrid      = 'off';
    hAx(i).YMinorGrid      = 'off';
end

saveas (hF,SaveName,'epsc')
disp([FileName,' saved to:'])
disp(PathName)

for i = 1:NAx
    hAx(i).XGrid           = 'on';
    hAx(i).YGrid           = 'on';
    hAx(i).XMinorGrid      = 'on';
    hAx(i).YMinorGrid      = 'on';
end