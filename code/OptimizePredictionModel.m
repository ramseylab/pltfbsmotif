function [xbest, ...
          score]=OptimizePredictionModel(trackData, ...
                                         bindingSites, ...
                                         tfNames, ...
                                         motifTrackInd, ...
                                         fprCutoffForROCArea, ...
                                         optimOptions, ...
                                         x0)

F = length(trackData);
T = size(trackData{1}, 2);
X = size(trackData{1}, 1);
if X ~= size(bindingSites, 1)
  error 'inconsistent vector sizes';
end
if F ~= size(bindingSites, 2)
  error 'inconsistent matrix size for num factors';
end
if F ~= length(tfNames)
  error 'inconsistent length of tfNames cell array';
end

baseWeight = 0.01;

if nargin < 7 || isempty(x0)
  x0=nan((3*T-1),1);
  fprintf(1, 'no starting parameter values provided\n\n');
else
  if length(x0) ~= 3*T-1
    error('incorrect size for x0');
    fprintf(1, 'x0 before filling in unknown parameters:\n');
    x0
  end
end
indsFill = find(isnan(x0));

% precompute the column-wise max of the trackData matrix, since
% these max values don't change during the optimization
trackDataMin = inf(1,T);
trackDataMax = -inf(1,T);
for f=1:F
  trackDataMin = min(min(trackData{f}), trackDataMin);
  trackDataMax = max(max(trackData{f}), trackDataMax);
end

if length(indsFill) > 0
  rangeCenters = zeros(T,1);
  rangeWidths = zeros(T,1);
  weights = baseWeight*ones(T,1);
  notMotifInds = setdiff(1:T, motifTrackInd);
  weights(motifTrackInd) = 1-sum(weights(notMotifInds));
  indepWeights0 = x0((2*T+1):(3*T-1));
  indsNotNan = find(~isnan(indepWeights0));
  weights(notMotifInds(indsNotNan)) = indepWeights0(indsNotNan);
  sumWeights = sum(abs(weights));
  if sumWeights > 1
    weights = weights/sumWeights;
  end
  indepWeights = weights(notMotifInds);
  
  rangeMaxs = trackDataMax;
  rangeMins = trackDataMin;
  
  rangeCenters(motifTrackInd) = 0.99;
  rangeWidths(motifTrackInd) = 0.01;
  tFill = indsFill(indsFill <= T);
  
  x0(indsFill(indsFill <= 2*T)) = [rangeMins(tFill); ...
                                  rangeMaxs(tFill)];
  x0((2*T+1):(3*T-1)) = indepWeights;
end

fprintf(1, 'x0 after filling in unknown parameters:\n');
x0

% precompute the list of unique binding site codes, since it
% doesn't change during the optimization
uniqueBindingSiteCodes = cell(1,F);
for f=1:F
  uniqueBindingSiteCodes{f} = unique(bindingSites(:,f));
end


%profile on; %:DEBUG:

P = 3*T - 1;

u = zeros(P,1);
v = zeros(P,1);

u(1:T) = trackDataMin;
v(1:T) = trackDataMax;
u((T+1):2*T) = trackDataMin;
v((T+1):2*T) = trackDataMax;
u((2*T+1):(3*T-1)) = -1;
v((2*T+1):(3*T-1)) = 1;


ncall = optimOptions(1);
fglob = optimOptions(2);
ftol = optimOptions(3);

initialScore = ComputePerformanceScore(x0', ...
                                       trackData, ...
                                       bindingSites, ...
                                       tfNames, ...
                                       fprCutoffForROCArea, ...
                                       [], ...
                                       motifTrackInd, ...
                                       uniqueBindingSiteCodes, ...
                                       trackDataMin, ...
                                       trackDataMax);
fprintf(1, 'initial score: %f\n', initialScore);

outputFileName = ['/tmp/' uuidgen() '_SNOBFIT.mat'];
xbest = ConstrainedOptSNOBFIT(@(x) ComputePerformanceScore(x, ...
                                              trackData, ...
                                              bindingSites, ...
                                              tfNames, ...
                                              fprCutoffForROCArea, ...
                                                  [], ...
                                              motifTrackInd, ...
                                              uniqueBindingSiteCodes, ...
                                              trackDataMin, ...
                                              trackDataMax), ...
                 u, ...
                 v, ...
                 x0, ...
                 round(ncall/4), ...
                 fglob, ...
                 ftol, ...
                 outputFileName);
delete(outputFileName);
xbest

optimOptions = optimset('Display','iter', ...
                        'MaxFunEvals', ncall, ...
                        'TolFun', ftol);

[xbest, score] = fminsearch(@(x) ComputePerformanceScore(x, ...
                                                  trackData, ...
                                                  bindingSites, ...
                                                  tfNames, ...
                                                  fprCutoffForROCArea, ...
                                                  [], ...
                                                  motifTrackInd, ...
                                                  uniqueBindingSiteCodes, ...
                                                  trackDataMin, ...
                                                  trackDataMax), ...
                   xbest, ...
                   optimOptions);

fprintf(1, 'Best score: %f\n', score);

end

