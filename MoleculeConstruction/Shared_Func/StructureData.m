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
   end
   
   properties
       Nmodes
       LocAlphaM
   end
   
   methods
       function nmodes = get.Nmodes(obj)
           nmodes = size(obj.LocFreq,1);
       end
       
      function locAlphaM = get.LocAlphaM(obj)
           locAlphaM = reshape(obj.LocAlpha,[],3,3);
       end
   end

end