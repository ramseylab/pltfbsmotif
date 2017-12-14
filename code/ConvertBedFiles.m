function ConvertBedFile(inputBedFile, ...
                        inputBarFile, ...
                        genomeRegionsBedFile, ...
                        outputBedFile, ...
                        featureSize)

[barChroms, ...
 barData]=readbar(inputBarFile);

[chromNames, ...
 startCoords, ...
 endCoords, ...
 featNames, ...
 peakHeights] = textread(inputBedFile, ...
                         '%s %d %d %s %f', ...
                         'headerlines', 1);

[chromNamesGR, ...
 startCoordsGR, ...
 endCoordsGR] = textread(genomeRegionsBedFile, ...
                         '%s %d %d', ...
                         'headerlines', 1);

fid=fopen(outputBedFile, 'w+');
F = length(chromNames);
for f=1:F
  chromName = chromNames{f};
  
  % get the start and end coordinates of the peak region
  startCoord = startCoords(f);
  endCoord = endCoords(f);
  
  % find the genome regions for the chromosome on which this peak is
  % located
  chromInds = find(strcmp(chromNamesGR, chromName));
  G = length(chromInds);
  
  % get the peak BAR file data for this chromosome
  barChromInd = find(strcmp(barChroms, chromName));
  if isempty(barChromInd)
    error 'unknown chromosome';
  end
  chromBarData = barData{barChromInd};
  chromBarCoords = chromBarData(:,1);
  chromBarInds = find(chromBarCoords >= startCoord-1 & ...
                      chromBarCoords < endCoord);
  
  % get the bar file sample at which the signal (within the peak region) is maximal
  [maxValue, maxValueInd] = max(chromBarData(chromBarInds,2));
  midCoord = chromBarCoords(chromBarInds(maxValueInd)) + 1;
  % =======================================================================
  % NOTE:  this was the original method, used for the Oct. 2009 analysis:
  %  startCoord = max(min(midCoord - 0.5*featureSize,  startCoord),1);
  % =======================================================================
  startCoord = max(midCoord - 0.5*featureSize, 1);
  endCoord = midCoord + 0.5*featureSize;
  for g=1:G
    startCoordGR = startCoordsGR(chromInds(g));
    endCoordGR = endCoordsGR(chromInds(g));
    if IntervalOverlap(startCoord, endCoord, startCoordGR, endCoordGR) > 0
      startCoordTrunc = max(startCoordGR, startCoord);
      endCoordTrunc = min(endCoordGR, endCoord);
      fprintf(fid, '%s\t%d\t%d\t%s\t%f\n', ...
              chromName, ...
              startCoordTrunc, ...
              endCoordTrunc, ...
              featNames{f}, ...
              peakHeights(f));
      break;
    end
  end
end
fclose(fid);
