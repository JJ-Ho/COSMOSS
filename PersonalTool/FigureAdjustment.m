PathName = '~/Desktop/TMP';
hF = figure(2);
% FigType = 'Structure';
% FigType = 'TS_IR';
FigType = 'TS_Raman';
% FigType = 'Intensity';
% FigType = 'O_Space'; Mode_Ind = '5';
% FigType = 'Refined_Space'; Mode_Ind = '1-5';
% FigType = 'Refind_Structure';
% FigType = '2DSFG';

Slide = 6;
PHI   = 0;
PSI   = 0;
THETA = 0;

Slide_Str = num2str(Slide);
PHI_Str   = num2str(PHI);
PSI_Str   = num2str(PSI);
THETA_Str = num2str(THETA);

FigTitleBase = ['PB-R6S3-X',Slide_Str];
FNameBase = ['PB_R6S3_X',Slide_Str];
FNameOrient = [FNameBase,'_',PHI_Str,'-',PSI_Str,'-',THETA_Str];
Orientaion = [', (\phi,\psi,\theta) = (',PHI_Str,',',PSI_Str,',',THETA_Str,')'];

hAx = findobj(hF,'Type','Axes');

%% Save picture
switch FigType
    case 'Structure'
        % Structure
        FigureName_Structure = [FNameOrient,'_Structure'];
        SaveFigures(hF,PathName,FigureName_Structure) 
        
    case 'TS_Raman'
        % Transition Strength Raman
        hAx = hAx(1);
        Orig_Title = hAx.Title.String;
        hAx.Title.String = [Orig_Title,Orientaion];

        FigureName_TS_Raman = [FNameBase,'_Raman'];

        SaveName = [PathName,'/',FigureName_TS_Raman];
        saveas (hF,SaveName,'png')

    case 'TS_IR'
        % Transition Strength IR
        hAx = hAx(1);
        Orig_Title = hAx.Title.String;
        hAx.Title.String = [Orig_Title,Orientaion];

        FigureName_TS_IR = [FNameBase,'_IR'];
        SaveName = [PathName,'/',FigureName_TS_IR];
        saveas (hF,SaveName,'png')
    
    case 'O_Space'
        % Orientation Space of a given mode
        Orig_Title = hAx.Title.String;
        hAx.Title.String = [FigTitleBase,', ',Orig_Title];

        FigureName_O_Space = [FNameBase,'_O-Space_Mode-',Mode_Ind];
        SaveFigures(hF,PathName,FigureName_O_Space)
        
    case 'Refined_Space'
        % Refined Orientation Space
        hAx.Title.String = [FigTitleBase,Orientaion,', Mode#',Mode_Ind];

        FigureName_RefinedSpace = [FNameOrient,'_RefinedSpace_Mode-',Mode_Ind];
        SaveFigures(hF,PathName,FigureName_RefinedSpace)
       
    case 'Refind_Structure'
        FigureName_RefinedStructure = [FNameOrient,'_Structure'];
        SaveFigures(hF,PathName,FigureName_RefinedStructure)
        
    case 'Intensity'
        FigureName_Int = [FNameBase,'_Int'];
        SaveFigures(hF,PathName,FigureName_Int)
        
    case '2DSFG'
        % 2DSFG
        hAx.XLim = [1500,1850];
        hAx.YLim = [1550,1800];
        Orig_Title = hAx.Title.String;
        hAx.Title.String = [Orig_Title,Orientaion];

        FigureName_2DSFG = [FNameOrient,'_2DSFG'];
        SaveFigures(hF,PathName,FigureName_2DSFG)
end