function AddBarFiles(inputFiles, ...
                     outputFile, ...
                     sampleInterval, ...
                     genomeChromSizes, ...
                     multiplier) 
%AddBarFiles - combines multiple BAR files
%
%  AddBarFiles(inputFiles, outputFile, 
%                            sampleInterval, genomeChromSizes,
%                            multiplier)
%
% Takes a cell array of BAR files, and adds them up at survey
% points along the chromosome, and multiplies the result by
% a multiplier. This is used for averaging the fragment counts
% from multiple fragment count files.  The resulting sample
% values are saved in a BAR file. NOTE: all input BAR files
% must be based on the same sampling interval (see 
% ConvertBedToBarFile.m).  Please note that this function
% is extremely memory-intensive, and should only be run
% on a 64-bit machine (at least when processing mammalian
% sequence data).
%
% Inputs:
%   inputFiles - a cell array containing the (full) names of
%    the BAR files to be combined; they must all have the
%   outputFile - the BAR file that will be written, containing
%    the combined (summed/rescaled) values of the input BAR files.
%   sampleInterval - the sampling interval of the intput BAR
%    files, in bp.  This needs to agree exactly with what was
%    used when the input BAR files were generated using
%    CnvtBedToBarFile.m.  In theory, an integer multiple of the
%    input BAR file sampling interval would also work.
%   genomeChromSizes - a struct array giving the size of each
%     chromosome, in bp.  The fields of the struct are:
%        * chromName:  the name of the chromosome (e.g., "chr2")
%        * chromSize:  the size of the chromosome, i bp
%     This data structure can be read in from a file, using the
%     function "GetGenomeChromSizes"
%   multiplier - a scalar by which the summed sample values
%     are rescaled before being saved as the output BAR file.
%     To do an average of the input files, just set this value
%     to the reciprocal of the number of input BAR files. 
%
% Outputs:
%
% Subfunctions:            none
% MAT-files required:      none
% See also:                none
%

% =================================================
% AddBarFiles.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================
G = length(genomeChromSizes);
chromNames = { genomeChromSizes(:).chromName };
chromSizes =  [ genomeChromSizes(:).chromSize ];

for g=1:G
  chromSize = chromSizes(g);
  chromCoords = (0:sampleInterval:(chromSize-1))';
  eval(sprintf('genomeChromData%d = zeros(length(chromCoords),1);', ...
               g));
end

I = length(inputFiles);
for i=1:I
% for each BAR file
  fileName= inputFiles{i};
  [seqNames, data] = readbar(fileName);
  C = length(seqNames);
  for c=1:C
    % for each chromosome
    seqName = seqNames{c};
    genomeChromInd = find(strcmp(chromNames, seqName));
    if isempty(genomeChromInd)
      error('unrecognized sequence name: %s in bar file: %s', seqName, ...
                    fileName);
    end
    seqData = data{c};
    chromSize = chromSizes(genomeChromInd);
    chromCoords = [0:sampleInterval:(chromSize-1)]';
    eval(sprintf('chromCounts = genomeChromData%d;', genomeChromInd));
    [intCoords, coordInds, seqInds] = intersect(chromCoords, seqData(:,1));
    chromCounts(coordInds) = chromCounts(coordInds) + ...
                             seqData(seqInds,2);
    eval(sprintf('genomeChromData%d = chromCounts;', genomeChromInd));
  end
end

for g=1:G
  eval(sprintf('genomeChromData%d = genomeChromData%d * multiplier;', ...
                               g, g));
end

outChromData = {};
for g=1:G
  chromSize = chromSizes(g);
  chromCoords = (0:sampleInterval:(chromSize-1))';
  eval(sprintf('chromCounts = genomeChromData%d;', g));
  chromInds = find(chromCounts);
  outChromData{g} = [chromCoords(chromInds) chromCounts(chromInds)];
end

CellToWriteBar(outputFile, ...
               chromNames, ...
               outChromData);


