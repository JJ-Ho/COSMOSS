function hF = SD_Draw(obj_SD,varargin)
    % Given the Structure data, this method will generate the molecular
    % structure with the proper figure plotting function
    hAx = 'New';
    if nargin > 1 
        hAx = varargin{:};
    end

    if isa(obj_SD.hPlotFunc,'function_handle')
        % check if triggered by a GUI interface
        if isempty(obj_SD.GUI_Inputs)
          obj_SD.GUI_Inputs.Debug = 'Debug';
          disp('There is no GUI_Inpus in the StructureData, so use default values for drawing...')
        end
        hF = obj_SD.hPlotFunc(hAx,obj_SD);
    else
        hF = '';
        disp('No @hPlotFunc defined, method "Draw" would not work...')
    end
end