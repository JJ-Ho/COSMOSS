classdef ColumnFormat 
   properties
      Name
      Format
      Width
      Data
      SortInd
   end
   
   methods
       
       function CF = ColumnFormat(Name,Format,Width)
           CF.Name    = Name;
           CF.Format  = Format;
           CF.Width   = Width;
           CF.Data    = [];
           CF.SortInd = [];
       end   
       
       function CF = ImportSortInd(CF,Data)
           % This method save the numerical data into Data and auto save
           % it to SortInd, too. When a column having string as Data, we
           % can update the "Data" structure later.
           CF.Data    = num2cell(Data); % trasform numeric data into strings
           CF.SortInd = abs(Data); % remove negative value for sorting
       end
           
       function MC = merge2cell(varargin)
           MC.Name    = {};
           MC.Format  = {};
           MC.Width   = {};
           MC.Data    = {}; % always as strings
           MC.SortInd = []; % must be numerical so use array instead of cell
           for i = 1:nargin
               MC.Name    = [MC.Name    ,varargin{i}.Name   ];
               MC.Format  = [MC.Format  ,varargin{i}.Format ];
               MC.Width   = [MC.Width   ,varargin{i}.Width  ];
               MC.Data    = [MC.Data    ,varargin{i}.Data   ];
               MC.SortInd = [MC.SortInd ,varargin{i}.SortInd];
               
               %if isnumeric(varargin{i}.Data)
               %    Next_Data = num2cell(varargin{i}.Data);
               %else
               %    Next_Data = varargin{i}.Data;
               %end     
           end
       end
       
   end
end