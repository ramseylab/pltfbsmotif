function MaskCopyBarFile(maskBedFile, trackBarFile, outputBarFile, gcs, binSize, operator, default)
% operator == 1:  insertion of value from trackBarFile into the
%                 corresponding "bucket" in maskBarFile; default 0
% operator == 2:  buckets in maskBarFile set to the distance to the
%                 nearest nonzero sample in trackBarFile


[chrsMask, startCoords, endCoords]=textread(maskBedFile, '%s %d %d', 'headerlines', 1);
% correct for UCSC "+1" end coordinates
endCoords = endCoords - 1;

chrsMaskUnique = unique(chrsMask);

[chrsTrack, dataTrack]=readbar(trackBarFile);

C = length(chrsMaskUnique);

gcsNames = { gcs.chromName };

outputChrs = cell(1,1);
outputData = cell(1,1);
chrCount = 1;
for c=1:C
  chrName = chrsMaskUnique{c};
  
  gcsInd = find(strcmp(gcsNames, chrName));
  if isempty(gcsInd)
    error 'could not find chromsome size';
  end

  chrSize = gcs(gcsInd).chromSize;

  chromBucketCoords = 0:binSize:(chrSize-1);
  B = length(chromBucketCoords);
  chromBucketMask = zeros(1,B);
  if 1==operator
    if nargin < 7
      chromBucketTrack = zeros(1,B);
    else
      chromBucketTrack = default * ones(1,B);
    end
  else
    chromBucketTrack = chrSize * ones(1,B);
  end
  chrsMaskInds = find(strcmp(chrName, chrsMask));
  I = length(chrsMaskInds);
  for i=1:I
    startInd = max(1, 1+floor((startCoords(chrsMaskInds(i))+0.5*binSize)/binSize));
    endInd = min(B, 1+floor((endCoords(chrsMaskInds(i))+0.5*binSize)/binSize));
    chromBucketMask(startInd:endInd) = 1;
  end
  
  trackChrInd = find(strcmp(chrsTrack, chrName));
  if ~isempty(trackChrInd)
    if length(trackChrInd) > 1
      error 'multiple chromosomes'
    end
    dataTrackChr = dataTrack{trackChrInd};
    P = size(dataTrackChr,1);
    for p=1:P
      switch operator
       case 1
        sampleInd = min(max(1, 1+floor((dataTrackChr(p,1)-1+0.5*binSize)/binSize)),B);
        chromBucketTrack(sampleInd) = dataTrackChr(p,2);
       case 2
        sampleCoord = dataTrackChr(p,1);
        dists = abs(chromBucketCoords-sampleCoord)/1000;
        chromBucketTrack = min(chromBucketTrack, dists);
       otherwise
        error 'invalid operator'
      end
    end
  end

  indsUse = find(chromBucketMask);
  U = length(indsUse);
  outputChrData = zeros(U,2);
  outputChrData(:,1) = chromBucketCoords(indsUse);
  outputChrData(:,2) = chromBucketTrack(indsUse);

  outputChrs{chrCount} = chrName;
  outputData{chrCount} = outputChrData;
  chrCount = chrCount + 1;
end

CellToWriteBar(outputBarFile, outputChrs, outputData);

