classdef StructureData
   properties
       XYZ
       AtomName
       COM
       
       LocCenter
       LocFreq
       LocAnharm
       LocMu
       LocAlpha
       
       StructModel
       FilesName
       
       Extra
   end
   
   properties
       LocAlphaM
       Nmodes
       NAtoms
       NStucture
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
   end

end