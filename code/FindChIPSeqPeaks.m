function numPeaks = FindChIPSeqPeaks(sampleBar, ...
                                     controlBar, ...
                                     absMin, ...
                                     outputFilePrefix, ...
                                     sampleName, ...
                                     peakDist, ...
                                     ratioMin, ...
                                     controlAbsMax)
% FindChIPSeqPeaks - analyze files of ChIP-seq counts to find peaks
%
%   numPeaks = FindChIPSeqPeaks(sampleBar, controlBar,
%                               absMin, outputFilePrefix, 
%                               sampleName, peakDist, 
%                               ratioMin, controlAbsMax)
%
% This program takes one or two BAR files containing counts of
% overlapping fragments at various positions throughout the genome,
% and identifies peak regions.  The counts may be real-valued,
% because they may have been normalized (which involves multiplying
% by a real-valued scale factor).  Peak regions are identified by
% finding positions where the "sample" BAR file height is above an
% absolute minimum threshold level, and optionally, above a
% specified multiplier times the fragment cound in the "control"
% BAR file at the same positions.  A peak region is defined as a
% set of positions where the above criteria are satisfied, of which
% adjacent positions are never more than a user-specified minimum
% distance apart (usually 200 bp).  The resulting peak list is
% ordered by chromosome and by coordinate, and saved as a both a
% BAR file and as a BED file (where the scores in the BED file are
% computed as 5 times the maximum fragment count for the peak
% region in the "sample" BAR file, up to a maximum score of 1000).
% The BAR files read by this program will typically have been
% generated using the "CnvtBedToBarFile.m" program.
%
% Inputs:
%
% sampleBar - a string variable containing the full name/location
%   of the BAR file of fragment counts from the target-specific IP
%   sample.  
% controlBar - a string variable containing the full name/location
%   of the BAR file of fragment counts from a negative control IP.
%   To disable the comparative criterion (sample vs. control), an
%   empty string can be passed as this argument; in that case, it
%   would be an error to supply a value for the "ratioMin" argument.
% absMin - a scalar (nonnegative integer) specifying the absolute
%   minimum fragent count (bar height) for the target-specific
%   sample, that must be observed at a given position, for that
%   position to be considered a "peak".
% outputFilePrefix - a string variable containing the full
%   name/location of an output file, to which a ".bed" or ".bar"
%   suffix is added as appropriate to generate the actual output
%   file name.  The filename portion of this string (i.e., the
%   portion after the last directory separator character, "/") is
%   used to generate the "name" field in the output BED file (if
%   there is no directory separator in the string, then the whole
%   string is used as the BED file "name" field.
% sampleName - a string identifier of the sample, to be embedded
%   in the output BED file track name.
% peakDist - a positive integer that specifies the minimum distance
%   between two points in the genome (at which the read-count is
%   above threshold), to be considered as belonging to two separate
%   ChIP-seq peaks.  The typical value used for this argument is
%   200.
% ratioMin - an optional floating-point variable (must be positive)
%   that specifies the minimum ratio of the number of overlapping
%   fragments in the target-specific IP sample, to the number of
%   overlapping fragments in the control sample.  Only specify this
%   argument if a control BAR file is also specified (with the
%   "controlBar" argument).  To omit this argument but still pass
%   the "controlAbsMax" argument, just pass an empty vector for
%   the "ratioMin" argument.
% controlAbsMax - an optional floating-point variable (must be
%   positive) that specifies the maximum value that a survey point
%   can have in the control sample, and still be considered a peak.
%
% Other m-files required:  readbar.m
%                          writebar.m
%                          IsAllowedChromosome.m
% Subfunctions:            none
% MAT-files required:      none
% See also:                CnvtBedToBarFile.m
%

% =================================================
% FindPeaksFromBarFiles.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================

if nargin > 6
  if isempty(ratioMin)
    ratioMin = 0;
  end
  if nargin > 7
    if isempty(controlAbsMax)
      controlAbsMax = Inf;
    end
  end
end

[sampleChrs, sampleData]=readbar(sampleBar);
if ~isempty(controlBar)
  if nargin < 7 
    warning 'ignoring control data, since no thresholds for using the control data were provided (see ratioMin and controlAbsMax arguments)';
  end
  
  [controlChrs, controlData]=readbar(controlBar);
  [intChrs, isam, ictr] = intersect(sampleChrs, controlChrs);
  sampleChrs = sampleChrs(isam);
  sampleData = sampleData(isam);
  controlData = controlData(ictr);
else
  if nargin > 6 && ratioMin > 0
    error 'need to provide a control BAR file when using the ratio test';
  end
end

C = length(sampleChrs);

chrNamesNumeric = cell(1,C);
for c=1:C
  chrName = sampleChrs{c};
  tokens = regexp(chrName, 'chr(\d+)(.*)$','tokens');
  if ~isempty(tokens)
    chromNum = tokens{1}{1};
    otherStuff = tokens{1}{2};
    chrName = sprintf('chr%0.4d%s', str2double(chromNum), otherStuff);
  end
  chrNamesNumeric{c} = chrName;
end

[sortedChrNamesNumeric,si] = sort(chrNamesNumeric);
sampleChrs = sampleChrs(si);
sampleData = sampleData(si);
if ~isempty(controlBar)
    controlData = controlData(si);
end

peakBarHeights = zeros(0,1);
peakBarCoords = zeros(0,1);
peakFileChromInds = zeros(0,1);

for c=1:C
  if IsAllowedChromosome(sampleChrs{c})
    sampleBarCoords = sampleData{c}(:,1);
    sampleBarHeights = sampleData{c}(:,2);
    indAbs = find(sampleBarHeights >= absMin);
    sampleBarCoords = sampleBarCoords(indAbs);
    sampleBarHeights = sampleBarHeights(indAbs);
    if ~isempty(controlBar) 
      controlBarCoords = controlData{c}(:,1);
      controlBarHeights = controlData{c}(:,2);
      [coordsInt, sampleInds, controlInds] = intersect(sampleBarCoords, ...
						       controlBarCoords);
      exclCoords = coordsInt(...
	  (sampleBarHeights(sampleInds) < ratioMin * controlBarHeights(controlInds)) | ...
	  (controlBarHeights(controlInds) > controlAbsMax));
    else
      exclCoords = []; 
    end
    [finalBarCoords, ia] = setxor(sampleBarCoords, intersect(sampleBarCoords, ...
						  exclCoords));
    finalBarHeights = sampleBarHeights(ia);
    if length(finalBarCoords) > 0
      peakBarHeights = [peakBarHeights; finalBarHeights];
      peakBarCoords = [peakBarCoords; finalBarCoords];
      peakFileChromInds = [peakFileChromInds; c*ones(length(finalBarCoords), ...
						     1)];
    end
  end
end

outputBarFile = [outputFilePrefix '.bar'];

writebar(outputBarFile, ...
	 sampleChrs(peakFileChromInds), ...
	 peakBarCoords, ...
	 peakBarHeights);

P = length(peakBarCoords);
outputBedFile = [outputFilePrefix '.bed'];
fid=fopen(outputBedFile, 'w+');
fprintf(fid, 'track name=%s description="ChIPSeq Data for %s" useScore=1\n', ...
	sampleName, ...
	sampleName);
numPeaks = 0;
if P > 0
  lastPeakStart = 1;
  lastPeakChrom = peakFileChromInds(1);
  for p=2:(P+1)
    if p==(P+1) || ...
         peakBarCoords(p) - peakBarCoords(p-1) > peakDist || ...
         peakFileChromInds(p) ~= lastPeakChrom 
      chromName = sampleChrs{lastPeakChrom};
      peakStartCoord = peakBarCoords(lastPeakStart);
      peakEndCoord = peakBarCoords(p-1);
      maxPeakHeight = max(peakBarHeights(lastPeakStart:(p-1)));
      numPeaks = numPeaks + 1;
      % need to add 1 to the "peakEndCoord" to comply with UCSC BED
      % file format
      fprintf(fid, '%s\t%d\t%d\t%s\t%f\n', ...
              chromName, ...
              peakStartCoord, ...
              peakEndCoord + 1, ...
              sampleName, ...
              min(1000, maxPeakHeight));
      if p < P+1
        lastPeakStart = p;
        lastPeakChrom = peakFileChromInds(p);
      end
    end
  end
end
fclose(fid);

