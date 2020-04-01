classdef DefaultFileChooser < FileChooser
    methods
        function [file,folder,status] = chooseFile(chooser,varargin)
            [file,folder,status] = uigetfile(varargin{:});
        end
    end
end