function OptimizeAndValidateModel(tracksToIncludeGeneric, ...
                                  tracksToIncludeFactorSpecific, ...
                                  outputMatFileName, ...
                                  initCondFileName, ...
                                  factorsToInclude, ...
                                  batchMode)

if nargin < 6
  batchMode = 0;
end

try
  
  tfMotifScanFileType = 'Combined';
  sequenceBackgroundType = 'B0';
  factors = {'P50';
             'P65';
             'ATF3';
             'IRF1';
             'CEBPd'};
  if nargin < 5 || isempty(factorsToInclude)
    factorsToInclude = ones(1,length(factors));
  end
  
  genericDataTrackFileName = 'TrackMatrix_Generic.mat';
  genomeChromSizesFileName = 'GenomeChromSizesM37.txt';
  bindingSiteDataFileName = 'TrackMatrix_BindingSites_50bp.mat';

  ncallBaseCoeff = 300;
  fglob = 0;
  ftol = 1e-4;
  fprCutoffForROCArea = 0.01;
  numFactorsToUseForOptimization = 4;
  numROCStepsForValidation = 100;
  
  % ============================================================================

  [pathstr, name]=fileparts(outputMatFileName);
  logFileName = [name '.log'];
  
  if ~batchMode
    diary(logFileName);
    fid=1;
  else
    fid=fopen(logFileName, 'w+');
  end
  
  if nargin > 3 && ~isempty(initCondFileName)
    fprintf(fid, 'x0 param file:      %s\n', initCondFileName);
    if strcmp('.txt',initCondFileName((end-3):end))
      x0=load(initCondFileName,'-ascii');
    else
      if strcmp('.mat', initCondFileName((end-3):end))
        initCondFileData = load(initCondFileName);
        x0 = initCondFileData.xbestAverage;
        if ~isempty(tracksToIncludeGeneric)
          warning('Disregarding the tracksToIncludeGeneric argument, using values from initial conditions file');
        end
        if ~isempty(tracksToIncludeFactorSpecific)
          warning('Disregarding the tracksToIncludeFactorSpecific argument, using values from initial conditions file');
        end
        tracksToIncludeGeneric = initCondFileData.tracksToIncludeGeneric;
        tracksToIncludeFactorSpecific = initCondFileData.tracksToIncludeFactorSpecific;
        clear 'initCondFileData';
      else
        error('unknown extension for initial condition file');
      end
    end
    if size(x0,1) == 1
      x0 = x0';
    end
  else
    x0=[];
    fprintf(fid, 'x0 param file:      (none)\n');
  end
  
  tracksToIncludeGenericInds = find(tracksToIncludeGeneric);
  tracksToIncludeFactorSpecificInds = find(tracksToIncludeFactorSpecific);
  factorsToIncludeInds = find(factorsToInclude);
  T = length(tracksToIncludeGenericInds) + ... 
      length(tracksToIncludeFactorSpecificInds);
  ncall = ceil(ncallBaseCoeff * T^2);

  optimOptions = [ncall, fglob, ftol];
  startTimeNum = now();
  startDateTime = datestr(startTimeNum);
  fprintf(fid, 'Processing started: %s\n', datestr(startTimeNum));

  currentWorkingDirectory = pwd();
  fprintf(fid, 'Working directory:  %s\n', currentWorkingDirectory);
  
  fprintf(fid, 'Output file:        %s\n', outputMatFileName);
  fprintf(fid, 'Max func calls:     %d\n', ncall);
  fprintf(fid, 'Function tolerance: %f\n', ftol);
  
  fprintf(fid, '\n');
  
  genomeChromSizes = GetGenomeChromSizes(genomeChromSizesFileName);
  motifTrackInd = length(find(tracksToIncludeGeneric)) + 1;
  fprintf(1, 'Loading generic track data file: %s\n', genericDataTrackFileName);
  genericData = load(genericDataTrackFileName);
  genericTrackData = genericData.trackData;
  if size(genericTrackData, 2) ~= length(tracksToIncludeGeneric)
    error 'invalid tracksToIncludeGeneric';
  end
  genericTrackFiles = genericData.trackFiles;
  
  genericChromInds = genericData.chromInds;
  chromNames = genericData.chromNames;
  chromCoords = genericData.chromCoords;
  clear 'genericData';

  bindingSiteData = load(bindingSiteDataFileName, ...
                         'trackData', ...
                         'trackFiles');
  bindingSites = bindingSiteData.trackData;
  bindingSites = bindingSites(:,factorsToIncludeInds);
  allTFNames = bindingSiteData.trackFiles;
  clear 'bindingSiteData';

  X = size(genericTrackData,1);

  F = length(factorsToIncludeInds);
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
    fprintf(fid, 'Loading data for factor: %s   File: %s\n', ...
            factorsUse{f}, factorSpecificTrackFileName);
    factorSpecificTrackFileNames{f} = factorSpecificTrackFileName;
    factorSpecificData = load(factorSpecificTrackFileNames{f});
    factorSpecificTrackDataFull = factorSpecificData.trackData;
    if f==1 && size(factorSpecificTrackDataFull, 2) ~= length(tracksToIncludeFactorSpecific)
      error 'invalid tracksToIncludeFactorSpecific';
    end
    factorSpecificTrackData{f} = factorSpecificTrackDataFull;
    factorSpecificTrackFileTrackNames{f} = factorSpecificData.trackFiles;
  end
  clear 'factorSpecificData';

  disp(' ');
  disp('Using generic tracks:');
  for i=1:length(tracksToIncludeGenericInds) 
    fprintf(1, '  %s\n', genericTrackFiles{tracksToIncludeGenericInds(i)});
  end
  disp(' ');
  disp('Using factor-specific tracks:');
  for i=1:length(tracksToIncludeFactorSpecificInds)
    fprintf(1, '  %s\n', factorSpecificTrackFileTrackNames{1}{tracksToIncludeFactorSpecificInds(i)});
  end
  fprintf(fid, '\n');
  
  optResults_rangeMins = zeros(N,T);
  optResults_rangeMaxs = zeros(N,T);
  optResults_weights = zeros(N,T);
  optResults_scoresOptimized = zeros(N,1);
  optResults_scoresValidation = nan(N,1);

  xbestSum = zeros(1,3*T - 1);

  for n=1:N

    indsToUseForOptimization = masterIndsToUseForOptimization(n,:);
    indsToUseForValidation = setdiff(1:F, indsToUseForOptimization);

    for f=1:numFactorsToUseForOptimization
      fUse = indsToUseForOptimization(f);
      factorSpecificTrackDataFull = factorSpecificTrackData{fUse};
      trackDataForOpt{f} = [genericTrackData(:, ...
                                             tracksToIncludeGenericInds) ...
                    factorSpecificTrackDataFull(:, ...
                                                tracksToIncludeFactorSpecificInds)];
      
    end

    [xbest, ...
     scoreOptimized] = OptimizePredictionModel(trackDataForOpt, ...
                                               bindingSites(:,...
                                                  indsToUseForOptimization), ...
                                               factorsUse(indsToUseForOptimization), ...
                                               motifTrackInd, ...
                                               fprCutoffForROCArea, ...
                                               optimOptions, ...
                                               x0);
    clear 'trackDataForOpt';
    
    optResults_rangeMins(n,:) = xbest(1:T);
    optResults_rangeMaxs(n,:) = xbest((T+1):2*T);
    indepWeights = xbest((2*T + 1):(3*T - 1));
    optResults_weights(n,setdiff(1:T, motifTrackInd)) = indepWeights;
    optResults_weights(n,motifTrackInd) = 1 - sum(abs(indepWeights));
    optResults_scoresOptimized(n) = scoreOptimized;
    
    if ~isempty(indsToUseForValidation)
      for f=1:numFactorsToUseForValidation
        fUse = indsToUseForValidation(f);
        factorSpecificTrackDataFull = factorSpecificTrackData{fUse};
        trackDataForVal{f} = [genericTrackData(:, ...
                                               tracksToIncludeGenericInds) ...
                    factorSpecificTrackDataFull(:, ...
                                                tracksToIncludeFactorSpecificInds)];
      end
      
      scoreValidation = ComputePerformanceScore(xbest, ...
                                                trackDataForVal, ...
                                                bindingSites(:,...
                                                  indsToUseForValidation), ...
                                                factorsUse(indsToUseForValidation), ...
                                                fprCutoffForROCArea, ...
                                                numROCStepsForValidation, ...
                                                motifTrackInd);
      
      optResults_scoresValidation(n) = scoreValidation;
      xbestSum = xbestSum + xbest/scoreValidation;
    end
    
    % flush the buffered "diary" output to disk (this step may only
    % be needed on glnxa64 architecture)
    if ~batchMode
      diary off;
      diary on;
    end
  end
  xbestAverage = xbestSum / sum(1./optResults_scoresValidation);

  endTimeNum = now();
  endDateTime = datestr(endTimeNum);
  elapsedTime = 24*60*(endTimeNum - startTimeNum);

  save(outputMatFileName, ...
       'optimOptions', ...
       'tfMotifScanFileType', ...
       'sequenceBackgroundType', ...
       'factors', ...
       'factorsToInclude', ...
       'tracksToIncludeGeneric', ...
       'tracksToIncludeFactorSpecific', ...
       'outputMatFileName', ...
       'genericDataTrackFileName', ...
       'genomeChromSizes', ...
       'factorSpecificTrackFileNames', ...
       'factorSpecificTrackFileTrackNames', ...
       'genericTrackFiles', ...
       'allTFNames', ...
       'bindingSiteDataFileName', ...
       'genomeChromSizesFileName', ...
       'currentWorkingDirectory', ...
       'startDateTime', ...
       'endDateTime', ...
       'elapsedTime', ...
       'optResults_weights', ...
       'optResults_rangeMins', ...
       'optResults_rangeMaxs', ...
       'optResults_scoresOptimized', ...
       'optResults_scoresValidation', ...
       'numFactorsToUseForOptimization', ...
       'xbestAverage', ...
       'ncallBaseCoeff', ...
       'fglob', ...
       'ftol', ...
       'motifTrackInd', ...
       'masterIndsToUseForOptimization', ...
       'numROCStepsForValidation');

catch
  % Errors need to be handled differently if the program is being
  % run in batch mode vs. interactively.  If in batch mode, print
  % the error message and exit MATLAB.  If interactive, just
  % rethrow the error and let MATLAB handle it normally.
  errStruct = lasterror;
  if nargin > 5 && batchMode
    errMessage = errStruct.message;
    errStack = errStruct.stack;
    L = length(errStack);
    for l=1:L
      fprintf(2, 'Error using ===> %s at %d\n%s\n', ...
	      errStack(l).name, ...
	      errStack(l).line, ...
	      errMessage);
    end
    exit;
  else
    rethrow(errStruct);
  end
end