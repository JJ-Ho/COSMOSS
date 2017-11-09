function ScoreGird = Score_Int(SigGrid,Target,Weight)

RelativeSig = bsxfun(@rdivide,SigGrid,SigGrid(:,1)); % Take ratio to the first signal
DistFromTarget = sqrt(bsxfun(@minus,RelativeSig,Target).^2);
WeightDist = bsxfun(@times,DistFromTarget,Weight);
ScoreGird  = 1./(sum(WeightDist,2) + 1);