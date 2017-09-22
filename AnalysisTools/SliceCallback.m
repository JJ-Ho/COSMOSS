function SliceCallback(hAx,hLs,Pos,h2D,Direction)
%% rename variables
hOsc  = hLs.hOsc;
X_Cut = Pos(1);
Y_Cut = Pos(3);

%% Find the X,Y,Z data
X = h2D.XData;
Y = h2D.YData;
Z = h2D.ZData;

%% update figure elements
N_format   = '%6.2f';

switch Direction
    case 'Y'
        [~,X_Cut_Ind] = min(abs(X - X_Cut));
        Osc_X = real(Z(:,X_Cut_Ind))';

        hOsc.YData  =  Osc_X;
        
        % Print XY values on the 2D figure title
        Origin_Str = hAx.Title.String;
        CVlaue_Str = regexprep(Origin_Str,'(-?[0-9]*\.[0-9]*', ['(',num2str(X_Cut,N_format)] );
        hAx.Title.String = CVlaue_Str;
        
    case 'X'     
        [~,Y_Cut_Ind] = min(abs(Y - Y_Cut));
        Osc_Y = real(Z(Y_Cut_Ind,:));
        
        hOsc.YData  =  Osc_Y;
        
        % Print XY values on the 2D figure title
        Origin_Str = hAx.Title.String;
        CVlaue_Str = regexprep(Origin_Str,'-?[0-9]*\.[0-9]*)', [num2str(Y_Cut,N_format),')'] );
        hAx.Title.String = CVlaue_Str;
end
