function MakeGenomeRegionsBedFile(startSiteTxtFile, ...
                                  genomeRegionFilePrefix, ...
                                  gcs, ...
                                  windowSize, ...
                                  binSize)

[chrs, coords]=textread(startSiteTxtFile, '%s %d');

gcs = GetGenomeChromSizes('../GenomeChromSizesM37.txt');
C = length(gcs);

data = cell(1,18);
chromNames = cell(1,18);
totalBins = 0;
chrCount = 1;
fid=fopen([genomeRegionFilePrefix '.bed'],'w+');
fprintf(fid, 'track name=GenomeRegions description=\"Genome Regions for SVM\"\n');
for c=1:C
  chrName = gcs(c).chromName
  chrLength = gcs(c).chromSize;
  bucketCoords = 0:binSize:(chrLength-1);
  B = length(bucketCoords);
  bucketVals = zeros(1,B);
  chrIndsSS = find(strcmp(chrs, chrName));
  S = length(chrIndsSS);
  for s=1:S
    coordSS = coords(chrIndsSS(s));
    minCoord = max(coordSS-windowSize,1);
    maxCoord = min(chrLength, coordSS+windowSize);
    minBucketInd = max(1+floor((minCoord+0.5*binSize)/binSize),1);
    maxBucketInd = min(B, 1+floor((maxCoord+0.5*binSize)/binSize));
    bucketVals(minBucketInd:maxBucketInd) = 1;
  end
  nzBucketInds = find(bucketVals);
  N = length(nzBucketInds);
  lastInd = 1;
  for n=2:N
    if nzBucketInds(n) > nzBucketInds(n-1)+1 || n==N
      if n < N
        endi = n-1;
      else
        endi = N;
      end
      fprintf(fid, '%s\t%d\t%d\n', chrName, ...
              bucketCoords(nzBucketInds(lastInd)), ...
              bucketCoords(nzBucketInds(endi))+1);
      lastInd = n; 
    end
  end
  if N > 0
    totalBins = totalBins + N;
    chrData = zeros(N,2);
    chrData(:,1) = bucketCoords(nzBucketInds);
    chrData(:,2) = 1;
    data{chrCount} = chrData;
    chromNames{chrCount} = chrName;
    chrCount = chrCount + 1;
  end
end
totalBins
fclose(fid);
CellToWriteBar([genomeRegionFilePrefix '.bar'], chromNames, data);


