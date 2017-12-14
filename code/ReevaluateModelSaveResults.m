function ReevaluateModelSaveResults(optResultsFile)

genericDataTrackFileName = '../../FinalTracks/TrackMatrix_Generic.mat';
bindingSiteDataFileName = '../../BindingSites/TrackMatrix_BindingSites.mat';
tfMotifScanFileType = 'Combined';
sequenceBackgroundType = 'B0';
fprCutoffForROCArea = 0.01;

optResultsFileData = load(optResultsFile);
xbest = optResultsFileData.xbestAverage;
N = length(xbest);
T = (N+1)/3;
tracksToIncludeGeneric = optResultsFileData.tracksToIncludeGeneric;
numFactorsToUseForOptimization = optResultsFileData.numFactorsToUseForOptimization;
motifTrackInd = length(find(tracksToIncludeGeneric)) + 1;

scoresValid = optResultsFileData.optResults_scoresValidation;
modelInds = find(scoresValid >= 5);

bindingSiteData = load(bindingSiteDataFileName, ...
                       'trackData', ...
                       'trackFiles');
bindingSites = bindingSiteData.trackData;
clear 'bindingSiteData';
factorsToInclude = optResultsFileData.factorsToInclude;
allFactorsToIncludeInds = find(factorsToInclude);
F = length(allFactorsToIncludeInds);
perms = nchoosek(allFactorsToIncludeInds,numFactorsToUseForOptimization);

genericTrackFileData = load(genericDataTrackFileName);
genericTrackData = genericTrackFileData.trackData;

perms

for modelInd = modelInds'
  
modelInd
  weights = optResultsFileData.optResults_weights(modelInd,:);
  indepWeights = zeros(1,T-1);
  indepWeights(:) = weights(setdiff(1:T,motifTrackInd));
  xbest = [optResultsFileData.optResults_rangeMins(modelInd,:) ...
           optResultsFileData.optResults_rangeMaxs(modelInd,:), ...
           indepWeights];

  factors = optResultsFileData.factors;
  tracksToIncludeFactorSpecific = optResultsFileData.tracksToIncludeFactorSpecific;
  tracksToIncludeGenericInds = find(tracksToIncludeGeneric);
  tracksToIncludeFactorSpecificInds = find(tracksToIncludeFactorSpecific);

  factorsToIncludeInds = setdiff(allFactorsToIncludeInds, perms(modelInd,:))

  factorsUse = factors(factorsToIncludeInds);
  trackData = cell(1,length(factorsToIncludeInds));
  for f=1:length(factorsToIncludeInds)
    factorSpecificTrackFileName = sprintf('../../FinalTracks/TrackMatrix_%s_%s_%s.mat', ...
                                          factorsUse{f}, ...
                                          tfMotifScanFileType, ...
                                          sequenceBackgroundType);
    factorSpecificTrackFileData = load(factorSpecificTrackFileName);
    factorSpecificTrackData = factorSpecificTrackFileData.trackData;
    trackData{f} = [genericTrackData(:,tracksToIncludeGenericInds) ...
                    factorSpecificTrackData(:,tracksToIncludeFactorSpecificInds)];
  end


  xbest

  score = ComputePerformanceScore(xbest, ...
                                  trackData, ...
                                  bindingSites(:,factorsToIncludeInds), ...
                                  factorsUse, ...
                                  fprCutoffForROCArea, ...
                                  motifTrackInd);

  optResultsFileData.optResults_scoresValidation(modelInd) = score;
end


outFile = ['Fixed/' optResultsFile];

save(outFile, '-struct', 'optResultsFileData');


