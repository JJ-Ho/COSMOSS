function obj_comb2 = SD_Comb2(obj1,obj2)

    XYZ_1       = obj1.XYZ;
    AtomName_1  = obj1.AtomName;
    LocCenter_1 = obj1.LocCenter;
    LocFreq_1   = obj1.LocFreq;
    LocAnharm_1 = obj1.LocAnharm;
    LocMu_1     = obj1.LocMu;
    LocAlpha_1  = obj1.LocAlpha;
    
    
    XYZ_2       = obj2.XYZ;
    AtomName_2  = obj2.AtomName;
    LocCenter_2 = obj2.LocCenter;
    LocFreq_2   = obj2.LocFreq;
    LocAnharm_2 = obj2.LocAnharm;
    LocMu_2     = obj2.LocMu;
    LocAlpha_2  = obj2.LocAlpha;
    

    XYZ        = [      XYZ_1;       XYZ_2];
    AtomName   = [ AtomName_1;  AtomName_2];
    LocCenter  = [LocCenter_1; LocCenter_2];
    LocFreq    = [  LocFreq_1;   LocFreq_2];
    LocAnharm  = [LocAnharm_1; LocAnharm_2];
    LocMu      = [    LocMu_1;     LocMu_2];
    LocAlpha   = [ LocAlpha_1;  LocAlpha_2];
    
    FilesName_1 = obj1.FilesName;
    FilesName_2 = obj2.FilesName;
    FilesName   = ['Comb2: ',FilesName_1,' & ',FilesName_2];
    
    % will be upgrade to real COM once I expand the AtomProperties.mat with the atomic mass
    COM_1 = obj1.COM;
    COM_2 = obj2.COM;
    COM = mean([COM_1;COM_2]); 
    
    % Need to merge the same filed of any?
    Extra_1 = obj1.Extra;
    Extra_2 = obj2.Extra;
    Extra = [Extra_1,Extra_2];

    obj_comb2 = StructureData;
    obj_comb2.StructModel = 5;
    obj_comb2.XYZ         = XYZ;
    obj_comb2.AtomName    = AtomName;
    obj_comb2.LocCenter   = LocCenter;
    obj_comb2.LocFreq     = LocFreq;
    obj_comb2.LocAnharm   = LocAnharm;
    obj_comb2.LocMu       = LocMu;
    obj_comb2.LocAlpha    = LocAlpha;
    
    obj_comb2.FilesName   = FilesName;
    obj_comb2.COM         = COM;
    obj_comb2.Extra       = Extra;
    
    
end