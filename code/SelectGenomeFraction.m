function indsToUseForOptimization=SelectGenomeFraction(genomeFracToUseForTraining, ...
                                                  X, ...
                                                  chromNames, ...
                                                  chromCoords, ...
                                                  genericChromInds, ...
                                                  genomeChromSizes, ...
                                                  bindingSites)

indsToUseForOptimization = zeros(1,ceil(0.7*genomeFracToUseForTraining*X));
numIndsAdded = 0;
C = length(chromNames);
for c=1:C
  chromSizeInd = find(strcmp({genomeChromSizes.chromName}, chromNames{c}));
  if isempty(chromSizeInd)
    error 'could not find chrom size';
  end
  chromSize = genomeChromSizes(chromSizeInd).chromSize;
  chromInds = find(genericChromInds == c);
  maxBindingValues = max(bindingSites(chromInds,:),[],2);
  indsNoBinding = find(~maxBindingValues);
  coordsNoBinding = chromCoords(chromInds(indsNoBinding));
  selectHighFractOfChrom = floor(2*rand());
  if selectHighFractOfChrom
    chromFrac = 1-genomeFracToUseForTraining;
  else
    chromFrac = genomeFracToUseForTraining;
  end
  indAbove = find(coordsNoBinding >= chromFrac * chromSize);
  [coordTrans, cti] = min(coordsNoBinding(indAbove));
  if ~isempty(coordTrans)
    indTrans = indsNoBinding(indAbove(cti));
    if selectHighFractOfChrom
      indsAdd = chromInds(indTrans:length(chromInds));
    else
      indsAdd = chromInds(1:indTrans);
    end
    indsToUseForOptimization((numIndsAdded+1):(numIndsAdded+length(indsAdd))) = ...
        indsAdd';
    numIndsAdded = numIndsAdded + length(indsAdd);
  end
end

