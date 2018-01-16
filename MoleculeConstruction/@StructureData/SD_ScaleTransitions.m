function obj_S = SD_ScaleTransitions(obj_SD,Scaling)
obj_S = SD_Copy(obj_SD);

obj_S.Scaled_LocMu    = Scaling.*obj_SD.LocMu;
obj_S.Scaled_LocAlpha = Scaling.*obj_SD.LocAlpha;