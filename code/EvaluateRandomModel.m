function EvaluateRandomModel

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

  fprCutoffForROCArea = 0.01;
  numFactorsToUseForOptimization = 4;


  % ============================================================================


  factorsToIncludeInds = 1:length(factors);


  x0=[];
  
 
  bindingSiteData = load(bindingSiteDataFileName, ...
                         'trackData', ...
                         'trackFiles');
  bindingSites = bindingSiteData.trackData;
  bindingSites = bindingSites(:,factorsToIncludeInds);
  allTFNames = bindingSiteData.trackFiles;
  clear 'bindingSiteData';

  F = size(bindingSites,2);
  optResults_scoresValidation = zeros(F,1);
  numSteps = 100;
  
  for f=1:F
    
    optResults_scoresValidation(f) = ComputePerformanceScore([], ...
                                                  [], ...
                                                  bindingSites(:,f), ...
                                                  factors(f), ...
                                                  fprCutoffForROCArea, ...
                                                  numSteps)
  end
  
  tracksToIncludeGeneric = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
  tracksToIncludeFactorSpecific = [0 0 0 0 0 0 0 0];
  save 'Output_RandomModel.mat' 'optResults_scoresValidation' 'tracksToIncludeGeneric' 'tracksToIncludeFactorSpecific';
  
    
      
