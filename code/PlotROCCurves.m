function [allFPRVals, allTPRVals]=PlotROCCurves2(outputMatFiles, legendEntries, fprCutoffForROCArea)

  tfMotifScanFileType = 'Combined';
  sequenceBackgroundType = 'B0';
  factors = {'P50';
             'P65';
             'ATF3';
             'IRF1';
             'CEBPd'};
  
  genericDataTrackFileName = 'TrackMatrix_Generic.mat';
  genomeChromSizesFileName = 'GenomeChromSizesM37.txt';

  bindingSiteDataFileName = 'TrackMatrix_BindingSites_50bp.mat';

  genericData = load(genericDataTrackFileName);
  genericTrackData = genericData.trackData;
  X = size(genericTrackData,1);
  genericTrackFiles = genericData.trackFiles;
  genericChromInds = genericData.chromInds;
  chromNames = genericData.chromNames;
  chromCoords = genericData.chromCoords;
  numSteps = 100;
  
  fprVals = 0:(fprCutoffForROCArea/numSteps):fprCutoffForROCArea;
  
  M = length(outputMatFiles);

  genomeChromSizes = GetGenomeChromSizes(genomeChromSizesFileName);

  bindingSiteData = load(bindingSiteDataFileName, ...
                         'trackData', ...
                         'trackFiles');
  bindingSites = bindingSiteData.trackData;
  allTFNames = bindingSiteData.trackFiles;

  colorsUse = [0 0 0;
              0 1 0;
              1 0 0;
              0 0 1;
              0 1 1;
              1 0 1;
              1 1 0;
              0.25 0.6 0.25;
              0.5 0.5 0.5];
  
  outputMatFileData = load(outputMatFiles{2});
  factorsToInclude = outputMatFileData.factorsToInclude;
  factorsToIncludeInds = find(factorsToInclude);
  F = length(factorsToIncludeInds);
  numFactorsToUseForOptimization = outputMatFileData.numFactorsToUseForOptimization;
  
  bindingSites = bindingSites(:,factorsToIncludeInds);
  masterIndsToUseForOptimization = nchoosek(1:F,numFactorsToUseForOptimization);
  N = size(masterIndsToUseForOptimization,1);
  numFactorsToUseForValidation = F-numFactorsToUseForOptimization;
  
  factorsUse = factors(factorsToIncludeInds);
  factorSpecificTrackData = cell(1,F);
  factorSpecificTrackFileNames = cell(1,F);
  factorSpecificTrackFileTrackNames = cell(1,F);
  
  trackDataForOpt = cell(1,numFactorsToUseForOptimization);
  trackDataForVal = cell(1,numFactorsToUseForValidation);
      
  for f=1:F
    factorSpecificTrackFileName = sprintf('TrackMatrix_%s_%s_%s.mat', ...
                                          factorsUse{f}, ...
                                          tfMotifScanFileType, ...
                                          sequenceBackgroundType);
    factorSpecificTrackFileNames{f} = factorSpecificTrackFileName;
    factorSpecificData = load(factorSpecificTrackFileNames{f});
    factorSpecificTrackDataFull = factorSpecificData.trackData;
    factorSpecificTrackData{f} = factorSpecificTrackDataFull;
    factorSpecificTrackFileTrackNames{f} = factorSpecificData.trackFiles;
  end      
  
  Q = length(fprVals);
  allFPRVals = zeros(M,Q);
  allTPRVals = zeros(M,Q);
  
  for m=1:M
    
    tprVals = zeros(size(fprVals));
    
    outputMatFile = outputMatFiles{m}
    sumWeights = 0;
    
    outputMatFileData = load(outputMatFile);
    tracksToIncludeFactorSpecific = outputMatFileData.tracksToIncludeFactorSpecific;
    tracksToIncludeGeneric = outputMatFileData.tracksToIncludeGeneric;
    tracksToIncludeGenericInds = find(tracksToIncludeGeneric);

    if size(genericTrackData, 2) ~= length(tracksToIncludeGeneric)
      error 'invalid tracksToIncludeGeneric';
    end
    tracksToIncludeFactorSpecificInds = find(tracksToIncludeFactorSpecific);
    T = length(tracksToIncludeGenericInds) + ... 
        length(tracksToIncludeFactorSpecificInds);
    
    if T > 0
      optResults_weights = outputMatFileData.optResults_weights;
      optResults_rangeMins = outputMatFileData.optResults_rangeMins;
      optResults_rangeMaxs = outputMatFileData.optResults_rangeMaxs;
      motifTrackInd = outputMatFileData.motifTrackInd;
    else
      optResults_weights = [];
      optResults_rangeMins = [];
      optResults_rangeMaxs = [];
      motifTrackInd = [];
    end

    for n=1:N
      % for n=1:N
      indsToUseForOptimization = masterIndsToUseForOptimization(n,:);
      indsToUseForValidation = setdiff(1:F, indsToUseForOptimization);

      if T > 0
        xbest = zeros(1,(3*T - 1));
        xbest(1:T) = optResults_rangeMins(n,:);
        xbest((T+1):2*T) = optResults_rangeMaxs(n,:);
        indsNotMotifs = setdiff(1:T, motifTrackInd)
        xbest((2*T + 1):(3*T - 1)) = optResults_weights(n,indsNotMotifs);
      else
        xbest = [];
      end
      
      for f=1:numFactorsToUseForValidation
        fUse = indsToUseForValidation(f);
        factorSpecificTrackDataFull = factorSpecificTrackData{fUse};
        trackDataForVal{f} = [genericTrackData(:, ...
                                               tracksToIncludeGenericInds) ...
                    factorSpecificTrackDataFull(:, ...
                                                tracksToIncludeFactorSpecificInds)];
      end
      
      xbest
      
      [scoreValidation, ...
       fprs, ...
       tprs] = ComputePerformanceScore(xbest, ...
                                       trackDataForVal, ...
                                       bindingSites(:,...
                                                    indsToUseForValidation), ...
                                       factorsUse(indsToUseForValidation), ...
                                       fprCutoffForROCArea, ...
                                       numSteps, ...
                                       motifTrackInd);
      
      auroc = fprCutoffForROCArea*(1-scoreValidation)
      weight = length(find(bindingSites(:,indsToUseForValidation)));
      sumWeights = sumWeights + weight;
      [fprsUniq, fprInds] = unique(fprs);
      tprVals = tprVals + interp1(fprsUniq, tprs(fprInds), fprVals);
    end
    
    allFPRVals(m,:) = fprVals;
    allTPRVals(m,:) = tprVals/N;
    plot(fprVals, tprVals/N, 'Color', colorsUse(m,:));
    hold on;
  end

  
lh = legend(legendEntries, 'Location', 'SouthEast');
set(lh, 'FontSize', 9);
xlabel('False positive error rate');
ylabel('Sensitivity (1 - false negative error rate)');
end

