function MaskBarFile(sampleBar, ...
                     controlBar, ...
                     outputBarFile, ...
                     controlAbsMax)
% MaskBarFile - mask a BAR file using thresholded values of a 2nd file
%
%   MaskBarFile(sampleBar, controlBar, outputBarFile,
%               controlAbsMax);
%
% This program reads in a BAR file and a "mask" BAR file, and sets
% any samples of the first BAR file, that happen to coincide with a
% position where the mask BAR file has a value exceeding a global
% threshold, to zero.  The masked output is written to a new BAR file.
%
% Inputs:
%
% sampleBar - a string variable containing the full name/location
%   of the BAR file of fragment counts from the target-specific IP
%   sample.  
% controlBar - the file to be used for masking the "sampleBar"
%   file.
% outputBarFiel - the name of the output (masked BAR) file to be written
% controlAbsMax - the threshold for masking; if the value in the
% control bar file exceeds this threshold at a position, the
%   corresponding sample in the "sampleBar" file at this position is
%   set to zero.
%
% Other m-files required:  readbar.m
%                          writebar.m
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

[sampleChrs, sampleData]=readbar(sampleBar);
[controlChrs, controlData]=readbar(controlBar);
[intChrs, isam, ictr] = intersect(sampleChrs, controlChrs);
sampleChrs = sampleChrs(isam);
sampleData = sampleData(isam);
controlData = controlData(ictr);

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
  sampleBarCoords = sampleData{c}(:,1);
  sampleBarHeights = sampleData{c}(:,2);
  indAbs = find(sampleBarHeights > 0);
  sampleBarCoords = sampleBarCoords(indAbs);
  sampleBarHeights = sampleBarHeights(indAbs);
  controlBarCoords = controlData{c}(:,1);
  controlBarHeights = controlData{c}(:,2);
  [coordsInt, sampleInds, controlInds] = intersect(sampleBarCoords, ...
                                                   controlBarCoords);
  exclCoords = coordsInt((controlBarHeights(controlInds) > controlAbsMax));
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

clear 'controlBarCoords';
clear 'controlBarHeights';
clear 'sampleBarCoords';
clear 'sampleBarHeights';
clear 'finalBarCoords';
clear 'finalBarHeights';

writebar(outputBarFile, ...
	 sampleChrs(peakFileChromInds), ...
	 peakBarCoords, ...
	 peakBarHeights);

