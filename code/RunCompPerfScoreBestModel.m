function RunCompPerfScoreBestModel(optResultsFile, modelInd)

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
weights = optResultsFileData.optResults_weights(modelInd,:);
indepWeights = zeros(1,T-1);
indepWeights(:) = weights(setdiff(1:T,motifTrackInd));
xbest = [optResultsFileData.optResults_rangeMins(modelInd,:) ...
         optResultsFileData.optResults_rangeMaxs(modelInd,:), ...
         indepWeights];
factors = optResultsFileData.factors;
factorsToInclude = optResultsFileData.factorsToInclude;
tracksToIncludeFactorSpecific = optResultsFileData.tracksToIncludeFactorSpecific;
tracksToIncludeGenericInds = find(tracksToIncludeGeneric);
tracksToIncludeFactorSpecificInds = find(tracksToIncludeFactorSpecific);


clear 'optResultsFileData';

F = length(factors);
perms = nchoosek(1:F,numFactorsToUseForOptimization);
factorsToIncludeInds = setdiff(1:F, perms(modelInd,:))
%factorsToIncludeInds = find(factorsToInclude)


genericTrackFileData = load(genericDataTrackFileName);
genericTrackData = genericTrackFileData.trackData;

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
clear 'genericTrackFileData';
clear 'factorSpecificTrackFileData';

bindingSiteData = load(bindingSiteDataFileName, ...
                       'trackData', ...
                       'trackFiles');
bindingSites = bindingSiteData.trackData;
clear 'bindingSiteData';

xbest

score = ComputePerformanceScore(xbest, ...
                                trackData, ...
                                bindingSites(:,factorsToIncludeInds), ...
                                factorsUse, ...
                                fprCutoffForROCArea, ...
                                motifTrackInd);

score