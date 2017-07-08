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
      
      obj_T     = SD_Trans(obj,V)
      obj_R     = SD_Rot(obj,Phi,Psi,Theta)
      obj_comb2 = SD_Comb2(obj1,obj2)
      Dihedral  = SD_PeptideDihedral(obj)
      AP        = SD_AtomicProperties(obj)
      
   end

end