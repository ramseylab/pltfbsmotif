function [chromNames, data] = ScanAllSequences(sequences, seqLabels, matrix, background, gcs, outputBarFile, binSize)

S = length(sequences);

C = length(gcs);
chromNames = { gcs.chromName };
data = cell(1,C);
for c=1:C
  chromSize = gcs(c).chromSize;
  chromCoords = [0:binSize:(chromSize-1)];
  B = length(chromCoords);
  chromData = zeros(B,2);
  chromData(:,1) = chromCoords;
  chromData(:,2) = nan;
  data{c} = chromData;
end

for s=1:S
  seqLabel = seqLabels{s};
  labelTokens = regexp(seqLabel,'[^=]+=([^\:]+)\:(\d+)-(\d+) ','tokens');
  if length(labelTokens) == 0
    error 'unable to parse sequence label';
  end  
  chromName = labelTokens{1}{1};
  chromInd = find(strcmp(chromNames, chromName));
  if isempty(chromInd)
    error 'unable to find chromosome size';
  end
  chromSize = gcs(chromInd).chromSize;
  startCoord = str2num(labelTokens{1}{2}) - 1;
  endCoord = str2num(labelTokens{1}{3}) - 1;
  
  chromData = data{chromInd};
  B = size(chromData,1);
  
  sequence = sequences{s};
  
  seqScores = ScanSequenceForMatrix(sequence, matrix, background);
  
  seqCoords = startCoord:endCoord;
  
  [maxScore, indMax]=max(seqScores);
  
  binInds = min(B, max(1, 1 + floor((seqCoords + 0.5*binSize)/binSize)));
   
  uniqueBinInds = unique(binInds);
  U = length(uniqueBinInds);
  for u=1:U
    binInd = uniqueBinInds(u);
    binScore = max(seqScores(find(binInds==binInd)));
    chromData(binInd,2) = max(chromData(binInd,2), binScore);
  end
  
  data{chromInd} = chromData;
end

indsUse = [];
for c=1:C
  chromIndsUse = find(~isnan(data{c}(:,2)));
  if ~isempty(chromIndsUse)
    indsUse = [indsUse c];
    data{c} = data{c}(chromIndsUse,:);
  end
end

chromNames = chromNames(indsUse);
data = data(indsUse);

if ~isempty(outputBarFile)
  CellToWriteBar(outputBarFile, chromNames, data);
end


