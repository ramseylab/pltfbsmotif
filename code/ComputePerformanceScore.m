function [score, fprs, tprs] = ComputePerformanceScore(x, ...
                                         trackData, ...
                                         bindingSites, ...
                                         tfNames, ...
                                         fprCutoffForROCArea, ...
                                                  numSteps, ...
                                         motifTrackInd, ...
                                         uniqueBindingSiteCodes, ...
                                         trackDataMin, ...
                                         trackDataMax)

score = 0;
numTrim = 0;

N = length(x);

X = size(bindingSites,1);
F = size(bindingSites,2);

T = (N+1)/3;
if N > 0 && T ~= round(T)
  error('invalid size for the parameter vector for the model');
end

numTrueBindingSites = sum(bindingSites>0, 1);
if ~any(numTrueBindingSites)
  error('no binding sites were detected; cannot compute score');
end

scoreCode = '';

if size(x,2)==1
  x=x';
end

if N > 0
  rangeMins = x(1:T);
  rangeMaxs = x((T+1):2*T);
  
  % assign penalty if L1 measure of weights exceeds the value 1
  indepWeights = x((2*T + 1):(3*T - 1));
  sumWeights = sum(abs(indepWeights));
  if sumWeights > 1
    score = score + 100 + 100*(sumWeights - 1);
    scoreCode = [scoreCode 'L'];
  end

  % penalize negative range widths
  rangeWidths = rangeMaxs - rangeMins;
  rangeWidthsNegInds = find(rangeWidths < 0);
  if length(rangeWidthsNegInds) > 0
    score = score + 100 + 100*sum(abs(rangeWidths(rangeWidthsNegInds)));
    scoreCode = [scoreCode 'N'];
  end
  
  weights = zeros(1,T);
  weights(setdiff([1:T],motifTrackInd)) = indepWeights;
  weights(motifTrackInd) = 1 - sumWeights;

  if nargin > 8 && ~isempty(trackDataMin)
    inds = find((rangeMins < trackDataMin));
    if length(inds) > 0
      score = score + 100 + 100*sum(trackDataMin(inds) - rangeMins(inds));
      scoreCode = [scoreCode 'M'];
    end
  end
  
  if nargin > 9 && ~isempty(trackDataMax)
    inds = find((rangeMaxs > trackDataMax));
    if length(inds) > 0
      score = score + 100 + 100*sum(rangeMaxs(inds) - trackDataMax(inds));
      scoreCode = [scoreCode 'X'];
    end
  end
  
  if score > 0
    return;
  end
end

scoreCodeROC = '';
scores = zeros(1,F);

for f=1:F
  if nargin > 7 && ~isempty(uniqueBindingSiteCodes)
    uniqueBindingSiteCodesFac = uniqueBindingSiteCodes{f};
  else
    uniqueBindingSiteCodesFac = [];
  end

  if N > 0
    ypred = ComputePredictionScore(trackData{f}, ...
                                   weights, ...
                                   rangeMins, ...
                                   rangeMaxs, ...
                                   motifTrackInd);
  else
    ypred = rand(X,1);
  end
  
  [specScore, ...
   areaROC, ...
   specScoreCode, ...
   fprs, ...
   tprs] = ComputeAreaUnderROCCurve( ypred, ...
                                     bindingSites(:,f), ...
                                     fprCutoffForROCArea, ...
                                     uniqueBindingSiteCodesFac, ...
                                     numSteps);

  if (~isempty(specScoreCode == 'P')) && (nargin > 9) && ~isempty(trackDataMax)
    specScore = specScore + sum(trackDataMax - rangeMaxs);
  end
  scoreCodeROC = union(scoreCodeROC, specScoreCode);
  
%  fprintf('For factor %s, score is: %f  scoreCode: %s  areaROC: %f\n', tfNames{f}, specScore, specScoreCode, areaROC);
  scores(f) = specScore;
%  score = score + numTrueBindingSites(f) * specScore;   
end
sortScores = sort(scores);
score = mean(sortScores((numTrim+1):(F-numTrim)));
scoreCode = [scoreCode scoreCodeROC];

%fprintf(1, 'Score returned to optimizer is: %f\n', score);

end

