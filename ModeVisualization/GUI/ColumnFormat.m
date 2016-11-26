classdef ColumnFormat 
   properties
      Name
      Format
      Width
      Data
   end
   methods
       function CF = ColumnFormat(Name,Format,Width)
           CF.Name   = Name;
           CF.Format = Format;
           CF.Width  = Width;
           CF.Data   = [];
       end   
       function O = merge2cell(varargin)
           O.Name   = {};
           O.Format = {};
           O.Width  = {};
           O.Data   = [];
           for i = 1:nargin
               O.Name   = [O.Name  ,varargin{i}.Name];
               O.Format = [O.Format,varargin{i}.Format];
               O.Width  = [O.Width ,varargin{i}.Width];
               O.Data   = [O.Data  ,varargin{i}.Data];
           end
       end
   end
end