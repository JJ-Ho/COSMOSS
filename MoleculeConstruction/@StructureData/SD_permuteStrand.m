function obj_Suffled = SD_permuteStrand(obj_SD)
% Given ideal betasheet StructureData, this function randomly permute the
% local mode frequencies/anharmonicties as if the strand were shuffled. The
% XYZ position of the strands stay the same.
% Copyright Jia-Jung Ho, 2013-2020

% check if the input obj_SD is consructed by the ideal betasheet model
if isequal(obj_SD.hPlotFunc,@Plot_Betasheet_AmideI)
    N_Modes   = obj_SD.Nmodes;
    N_Residue = obj_SD.Extra.N_Residue;
    N_Strand  = obj_SD.Extra.N_Strand;
    L_Index   = obj_SD.Extra.L_Index;
    LocFreq   = obj_SD.LocFreq;
    LocAnharm = obj_SD.LocAnharm;
    
    p = randperm(N_Strand);
    LocFreq   = reshape(  LocFreq,N_Residue,N_Strand,[]);
    LocAnharm = reshape(LocAnharm,N_Residue,N_Strand,[]);
    
    LocFreq   = reshape(  LocFreq(:,p,:),N_Modes,[]);
    LocAnharm = reshape(LocAnharm(:,p,:),N_Modes,[]);
    
    L_Index   = reshape(L_Index,N_Residue,N_Strand);
    L_Index   = reshape(L_Index(:,p),N_Modes,[]);
    
    obj_Suffled                   = SD_Copy(obj_SD);
    obj_Suffled.LocFreq           = LocFreq;
    obj_Suffled.LocAnharm         = LocAnharm;
    obj_Suffled.Extra.permutIndex = p;
    obj_Suffled.Extra.L_Index     = L_Index;

else
    error('The SD_permuteStrand method only works for idealbeta sheet model...')
end