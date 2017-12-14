function BedToBarFileLabeledRegions(bedFile, genomeRegionsBarFile, binSize, outputBarFile)

[chromNames, ...
 startCoords, ...
 endCoords, ...
 featNames, ...
 peakHeights] = textread(bedFile, '%s %d %d %s %f');

startCoords = startCoords - 1;
endCoords = endCoords - 2;

[barChroms, barData]=readbar(genomeRegionsBarFile);
C = length(barChroms);
binCtr = 0;
for c=1:C
  chromData = barData{c};
  chromData(:,2) = 0;
  binCtr = binCtr + size(chromData,1);
  barData{c} = chromData;
end

v = zeros(binCtr, 1);

P = length(chromNames);
for p=1:P
  chromName = chromNames{p};
  chromInd = find(strcmp(barChroms, chromName));
  if isempty(chromInd)
    error 'could not find chromosome';
  end
  chromData = barData{chromInd};
  startCoord = startCoords(p);
  endCoord = endCoords(p);
  binInds = find(chromData(:,1) + 0.5*binSize >= startCoord & ...
                 chromData(:,1) - 0.5*binSize <= endCoord);
  chromData(binInds,2) = p;
  barData{chromInd} = chromData;
end

CellToWriteBar(outputBarFile, barChroms, barData);




