classdef StructureData < handle
   properties
       % Bare minimum properties that is needed for all the SD method 
       % start to work
       XYZ      % size = [NAtoms,3]
       AtomName % size = {NAtoms,1} 
   end
   
   properties
       LocCenter       % size = [Nmodes,3]
       LocFreq         % size = [Nmodes,1]
       LocAnharm       % size = [Nmodes,1]
       LocMu           % size = [Nmodes,3]
       LocAlpha        % size = [Nmodes,9]
       DiagDisorder    % size = [Nmodes,1]
       OffDiagDisorder % size = [Nmodes,1]
       
       FilesName
       
       Extra
       
       StructModel % this will be remove later
       Children    % this property is only used by Comb2, maybe redundent
       
       hPlotFunc 
       hGUIs         % this is used when calling the plot function, maybe redundent?
       hParseGUIFunc % this is used when calling the plot function, maybe redundent?
   end
   
   properties
       % These properties will be automatically update when the dependent
       % peoperty is assigned. But is is free to be change later.
       Scaled_LocMu    % for comb2 concentration scaling that only applys on the MuAlphaGen
       Scaled_LocAlpha % for comb2 concentration scaling that only applys on the MuAlphaGen
   end
   
   properties
       % These properties will be calculated when quaried
       Nmodes
       NAtoms
       NStucture
       LocAlphaM
       CoM
   end
   
   methods
      function Nmodes = get.Nmodes(obj)
           Nmodes = size(obj.LocFreq,1);
      end
      function NAtoms = get.NAtoms(obj)
           NAtoms = size(obj.XYZ,1);
      end
      function NStucture = get.NStucture(obj)
           NStucture = size(obj.XYZ,3);
      end       
      function locAlphaM = get.LocAlphaM(obj)
           % Vector version
           % [Aixx Aixy Aixz Aiyx Aiyy Aiyz Aizx Aizy Aizz] , size =[Mode_Num X 9]
           % 
           % Matrix representation
           % D3  --> D2
           % ^ [ Aixx Aixy Aixz ] 
           % | [ Aiyx Aiyy Aiyz ]
           % | [ Aizx Aizy Aizz ]
           locAlphaM = reshape(obj.LocAlpha,[],3,3);
      end      
      function CoM = get.CoM(obj)
           AP   = SD_AtomicProperties(obj);
           Mass = AP.Mass;
           CoM  = sum(bsxfun(@times,obj.XYZ,Mass),1)./sum(Mass);
      end
      function hF = Draw(obj,varargin)
          % Simplify the structure drawing syntex
          hAx = 'New';
          if nargin > 1 
              hAx = varargin{:};
          end
          if isa(obj.hPlotFunc,'function_handle')
              % check if triggered by a GUI interface
              if isempty(obj.hGUIs)
                  GUI_Input.Debug = 'Debug';
              else
                  GUI_Input = obj.hParseGUIFunc(obj.hGUIs);
              end
              hF = obj.hPlotFunc(hAx,obj,GUI_Input);
          else
              hF = '';
              disp('No @hPlotFunc defined, method "Draw" would not work...')
          end
      end
      
      % Automatically copy the non-scaled value to the scaled properties
      % when the value fisrt being assigned
      function set.LocMu(obj,Value)
          obj.LocMu = Value;
          obj.Scaled_LocMu = Value;
      end
      function set.LocAlpha(obj,Value)
          obj.LocAlpha = Value;
          obj.Scaled_LocAlpha = Value;
      end
      
      % Othe methods defined as a separated function
      AP         = SD_AtomicProperties(obj)
      obj_T      = SD_Trans(obj,V)
      obj_R      = SD_Rot(obj,R)
      obj_comb2  = SD_Comb2(obj1,obj2)
      Dihedral   = SD_PeptideDihedral(obj)
      obj_S      = SD_ScaleTransitions(obj,Scaling)
      obj_New    = SD_Copy(obj)
      obj_TN     = SD_TransN(obj,V,N)
      obj_Framed = SD_SetFrame(obj,Center_Ind,Z_Ind,XZ_Ind)
      Obj_AmideI = SD_GetAmideI(obj)
   end

end