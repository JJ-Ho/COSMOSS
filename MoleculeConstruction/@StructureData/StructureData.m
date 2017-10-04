classdef StructureData < handle
   properties
       XYZ
       AtomName
       
       LocCenter
       LocFreq
       LocAnharm
       LocMu
       LocAlpha
       
       StructModel
       FilesName
       
       Extra
       
       Children
       hPlotFunc
       hParseGUIFunc
       hGUIs
   end
   
   properties
       LocAlphaM
       Nmodes
       NAtoms
       NStucture
       CoM
   end
   
   methods
      function nmodes = get.Nmodes(obj)
           nmodes = size(obj.LocFreq,1);
      end
      function nAtoms = get.NAtoms(obj)
           nAtoms = size(obj.XYZ,1);
      end
      function nStucture = get.NStucture(obj)
           nStucture = size(obj.XYZ,3);
      end       
      function locAlphaM = get.LocAlphaM(obj)
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
              hF = obj.hPlotFunc(hAx,obj,obj.hParseGUIFunc(obj.hGUIs));
          else
              hF = '';
              disp('No @hPlotFunc defined, method "Draw" would not work...')
          end
      end
      
      AP        = SD_AtomicProperties(obj)
      obj_T     = SD_Trans(obj,V)
      obj_R     = SD_Rot(obj,Phi,Psi,Theta)
      obj_comb2 = SD_Comb2(obj1,obj2)
      Dihedral  = SD_PeptideDihedral(obj)
      obj_S     = SD_ScaleTransitions(obj,Scaling)
      obj_New   = SD_Copy(obj)
   end

end