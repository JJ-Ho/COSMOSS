function obj_S = SD_ScaleTransitions(obj,Scaling)
obj_S = SD_Copy(obj);

obj_S.Scaled_LocMu    = Scaling.*obj.LocMu;
obj_S.Scaled_LocAlpha = Scaling.*obj.LocAlpha;