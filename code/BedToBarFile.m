function BedToBarFile(bedFile, barFile, binSize, gcs, includeScores)
if nargin < 5
  includeScores=1;
end
if includeScores
  [chrs, startCoords, endCoords, names, scores]=textread(bedFile, ...
                                                  '%s %d %d %s %d');
else
  [chrs, startCoords, endCoords, names]=textread(bedFile, ...
                                                 '%s %d %d %s');
end  
% coorect for the "+1" end coordinates in the UCSC BED file
endCoords = endCoords - 1;
F = length(chrs);
uniqueChrs = unique(chrs);
C = length(uniqueChrs);
chrNamesGCS = { gcs.chromName };
outData = cell(1,C);
for c=1:C
  chrName = uniqueChrs{c};
  chrInd = find(strcmp(chrNamesGCS, chrName));
  if isempty(chrInd)
    error 'could not find chromosome';
  end
  chrSize = gcs(chrInd).chromSize;
  bucketCoords = 0:binSize:(chrSize-1);
  B = length(bucketCoords);
  bucketValues = zeros(1,B);
  chrFeatInds = find(strcmp(chrs, chrName));
  I = length(chrFeatInds);
  for i=1:I
    chrFeatInd = chrFeatInds(i);
    startCoord = startCoords(chrFeatInd);
    endCoord = endCoords(chrFeatInd);
    if includeScores
      score = scores(chrFeatInd);
    else
      score = 1;
    end
    startBucketInd = max(1, 1 + floor((startCoord+0.5*binSize)/binSize));
    endBucketInd = min(B, 1 + floor((endCoord+0.5*binSize)/binSize));
    indsRange = startBucketInd:endBucketInd;
    indsBelow = find(bucketValues(indsRange) < score);
    bucketValues(indsRange(indsBelow)) = score;
  end
  nzInds = find(bucketValues);
  N = length(nzInds);
  chrData = zeros(N,2);
  chrData(:,1) = bucketCoords(nzInds);
  chrData(:,2) = bucketValues(nzInds);
  outData{c} = chrData;
end

CellToWriteBar(barFile, uniqueChrs, outData);

