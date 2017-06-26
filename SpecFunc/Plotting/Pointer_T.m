function [] = Pointer_T(S)
% Put the cursor values at the second line of title

% pPrep variabes
hF  = S.hF;
hAx = S.hAx;

% Find old Pointer ui elements and delete if exist
hOld_Number = findobj(hF,'Tag','Pointer_Number');
if ishandle(hOld_Number)
    delete(hOld_Number)
end

% Find and replace the original title string with a new line for the cursor values
T_Str   = hAx.Title.String;
S.T_Str = T_Str;
T_C_Str = {T_Str,''};
hAx.Title.String = T_C_Str;
hAx.Title.FontSize = 12;

set(hF,'windowbuttonmotionfcn',{@fh_wbmfcn,S}) % Set the motion detector.

function [] = fh_wbmfcn(varargin)
% WindowButtonMotionFcn for the figure.
S = varargin{3};  % Get the structure.
PT = get(S.hAx,'CurrentPoint');  % The current point w.r.t the axes.
F = PT(1,1:2);

N_format   = '%6.2f';
CVlaue_Str = ['Current Point: (',num2str(F(1),N_format),', ',num2str(F(2),N_format),')'];
T_C_Str    = {S.T_Str,CVlaue_Str};

S.hAx.Title.String = T_C_Str;
