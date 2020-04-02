classdef DefaultFileChooser < FileChooser
    methods
        function [file,folder,status] = chooseFile(chooser,varargin)
            [file,folder,status] = uigetfile(varargin{:});
        end
    end
end

% PWD = pwd;
% PDB_Path = [PWD, '/StructureFiles/PDB/'];
% 
% [FilesName,PathName,~] = uigetfile({'*.pdb','PDB file'; ...
%                                     '*,*','All Files'},...
%                                     'MultiSelect','on',...
%                                     'Select inputs',PDB_Path);