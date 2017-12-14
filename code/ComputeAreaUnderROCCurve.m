function [score, ...
          rocCurveArea, ...
          scoreCode, ...
          fprs, ...
          tprs]=ComputeAreaUnderROCCurve(ypred, ...
                                         bindingSites, ...
                                         fprCutoffForROCArea, ...
                                         uniqueBindingSiteCodes, ...
                                         numSteps)

X = length(ypred);
if X ~= size(bindingSites, 1)
  error 'bindingSites size inconsistent with ypred';
end

scoreCode = '';
fprs = [];
tprs = [];

if nargin < 5 || isempty(numSteps)
  numSteps = 20;
end
minSamplesForROC = 3;

if nargin < 4 || isempty(uniqueBindingSiteCodes)
  uniqueBindingSiteCodes = unique(bindingSites);
end
uniqueBindingSiteCodes = uniqueBindingSiteCodes(uniqueBindingSiteCodes > 0);
B = length(uniqueBindingSiteCodes);

bindingSiteSize = length(find(bindingSites))/B;
tnEstimate = max(1, X - B*bindingSiteSize);

bindIndGlobal = find(bindingSites);
nonzeroBindingSites = bindingSites(bindIndGlobal);

ypredNZInds = find(ypred > 0);
numNZYpreds = length(ypredNZInds);
compQuantRangeMax = fprCutoffForROCArea;
compQuantRangeMin = 1/X;
if compQuantRangeMin >= compQuantRangeMax
  error('fprCutoffForROCArea is too small');
end

logCompQuantRangeMax = log(compQuantRangeMax);
logCompQuantRangeMin = log(compQuantRangeMin);

expStep = (logCompQuantRangeMax - logCompQuantRangeMin)/numSteps;
quantilesTry = 1-exp([logCompQuantRangeMin:expStep:logCompQuantRangeMax]);
if length(find(isnan(ypred))) > 0
  error 'found nans in ypred';
end
predValueQuants = quantile(ypred(ypredNZInds), quantilesTry);
if size(predValueQuants,1) ~= 1
  predValueQuants = predValueQuants';
end
predValueQuants=sort(unique(predValueQuants),1,'descend');
if length(predValueQuants) < minSamplesForROC
  % we couldn't get enough unique prediction values above an
  % estimated FPR cutoff of 0.01; need more diverse predictions
  score = 2*(minSamplesForROC-length(predValueQuants));
  scoreCode = [scoreCode 'P'];
  rocCurveArea=nan;
  return;
end

P = length(predValueQuants);
allPred = bsxfun(@ge, ypred, predValueQuants);
fns = zeros(1,P);
fps = zeros(1,P);
for b=uniqueBindingSiteCodes';
  fns = fns + ~any(allPred(bindIndGlobal(nonzeroBindingSites == b),:),1);
end
fps = sum(bsxfun(@and, allPred, ~bindingSites),1);
fprs = min(fps / tnEstimate, 1.0);
fnrs = min(fns / B, 1.0);
[fprs, fprsI] = sort(fprs);
fnrs = fnrs(fprsI);
tprs = 1 - fnrs;
indFprCutoff = max(find(fprs < fprCutoffForROCArea));
if length(indFprCutoff) == 0
  rocCurveArea = nan;
  score = 10*min(fprs)/fprCutoffForROCArea;
  scoreCode = [scoreCode 'F'];
  if isempty(score)
    error('empty score');
  end
  return;
end

indNext = min(find(fprs >= fprCutoffForROCArea));
if length(indNext > 0) 
  tprMax = tprs(indNext);
  fprMax = fprs(indNext);
else
  tprMax = 1;
  fprMax = 1;
end
fprs = [0 fprs(1:indFprCutoff) fprCutoffForROCArea];
xvals = [fprs(indFprCutoff) fprMax];
yvals = [tprs(indFprCutoff) tprMax];
  
tprInterp = interp1(xvals, ...
                    yvals, fprCutoffForROCArea);
tprs = [0 tprs(1:indFprCutoff) tprInterp];


rocCurveArea = trapz(fprs, tprs);
if rocCurveArea > fprCutoffForROCArea
  error 'roc curve area is too big';
end
score = max(fprCutoffForROCArea - rocCurveArea, 0.0)/fprCutoffForROCArea;
if isempty(score)
  error('empty score');
end
%fprintf(1, 'Area under ROC curve is: %f  num samples: %d\n', rocCurveArea, length(predValueQuants));

%fprintf(1, 'The score based on the ROC curve is: %f\n', score);
% minMaxDiff = rangeMins-rangeMaxs;
% minMaxDiffInd = find(minMaxDiff > 0);
% if length(minMaxDiffInd) > 0
%   score = score + sum(1./(minMaxDiff(minMaxDiffInd) + 0.001));
% %  fprintf(1, 'A penalty has been added; the score is now: %f\n', score);
% end

%figure;
%tvs = [0:(fprCutoffForROCArea/100):fprCutoffForROCArea];
%plot(fprs, tprs, '-b', ...
%     tvs, tvs, '--g');
%xlim([0 fprCutoffForROCArea]);
%ylim([0 0.25]);




