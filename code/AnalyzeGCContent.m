function AnalyzeGCContent(seqs, seqLabels, gcs, outBarFile, binSize)

%[seqs,seqLabels]=readfastaseqs('GenomeRegions_NotMasked.fa');
%gcs = GetGenomeChromSizes('../GenomeChromSizesM37.txt');
%binSize = 100;

C = length(gcs);
buckets = cell(1,C);
coords = cell(1,C);
chroms = cell(1,C);
for c=1:C
  chromName = gcs(c).chromName;
  chroms{c} = chromName;
  chromSize = gcs(c).chromSize;
  chromCoords = 0:binSize:(chromSize-1);
  coords{c} = chromCoords;
  B = length(chromCoords);
  buckets{c} = -1 * ones(1,B);
end

S = length(seqs);
for s=1:S
  seqLabel = seqLabels{s};
  labelTokens = regexp(seqLabel,'[^=]+=([^\:]+)\:(\d+)-(\d+) ','tokens');
  if length(labelTokens) > 0
    chrom = labelTokens{1}{1};
    chromInd = find(strcmp(chroms, chrom));
    if isempty(chromInd)
      error 'could not find data for chromosome';
    end
    chromCoords = coords{chromInd};
    B = length(chromCoords);
    startCoord = str2num(labelTokens{1}{2}) - 1;
    endCoord = str2num(labelTokens{1}{3}) - 1;
    seqCoords = startCoord:endCoord;
    binInds = min(max(1,1+floor((seqCoords + 0.5*binSize)/binSize)),B);
    uniqueBinInds = unique(binInds);
    U = length(uniqueBinInds);
    seq = seqs{s}; 
    if length(seq) ~= endCoord-startCoord+1
      error 'sequence length inconsistent with coordinate range';
    end
    chromBuckets = buckets{chromInd};
    for u=1:U
      binInd = uniqueBinInds(u);
      indsUse = find(binInds == binInd);
      seqUse = seq(indsUse);
      totalSeq = length(find(seqUse ~= 'N'));
      gcSeq = length(find(seqUse == 'G' | seqUse == 'C'));
      gcContent = gcSeq/totalSeq;
      chromBuckets(binInd) = gcContent;
    end
    buckets{chromInd} = chromBuckets;
  end
end

data = cell(1,C);
chromsUse = [];
for c=1:C
  chromBuckets = buckets{c};
  chromCoords = coords{c};
  nzBucketInds = find(chromBuckets >= 0);
  N = length(nzBucketInds);
  if N > 0
    chromData = zeros(N,2);
    chromData(:,1) = chromCoords(nzBucketInds);
    chromData(:,2) = chromBuckets(nzBucketInds);
    chromsUse = [chromsUse c];
    data{c} = chromData;
  end
end

CellToWriteBar(outBarFile, chroms(chromsUse), data(chromsUse));
