function AvgFragCountFiles(fragCountFiles, ...
                           fragFilesDir, ...
                           outputFile, ...
                           sampleInterval, ...
                           genomeChromSizes)
% AvgFragCountFiles - Averages multiple ChIP-seq fragment count files
%
% AvgFragCountFiles(fragCountFiles, fragFilesDir, outputFile, 
%                   sampleInterval, genomeChromSizes)
%
% This function reads in a set of Affymetrix BAR files that give
% the number of overlapping fragments at various positions in the
% genome, and computes the total average number of fragments at
% survey points established at fixed intervals along the
% chromosome.  The resulting average fragment counts are scaled to
% "tags per million" by dividing by the total number of fragments
% (across all fragment count files being averaged) and then
% multiplying by 1 million.
%
% Inputs:
%  fragCountFiles - a cell array containing the names of the BAR
%    files whose data points are to be averaged
%  fragFilesDir - a string containing the name of the directory
%    containing the extended fragment files, in BED format (the
%    total number of fragments in each BAR file is obtained from 
%    the line count of the corresponding BED file).
%  outputFile - the name of the output BAR file to be written
%  sampleInterval - the sampling interval (in bp) for the BAR file to
%    be written, representing the average of the input BAR files
%  gcs - a struct array containing the names and sizes (in bp) of
%    all of the chromosomes (see GetGenomeChromSizes.m)
%  
%  Uses:
%    AddBarFiles.m
%

% =================================================
% AvgFragCountFiles.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================

F = length(fragCountFiles);

numFragments = 0;
for f=1:F
  fileName = fragCountFiles{f};
  [pathstr, name]=fileparts(fileName);
  bedFileName=[fragFilesDir '/' name '.bed'];
  [status, result]=system(['cat ' bedFileName ' | wc -l | cut -f1']);
  numFragments = numFragments + str2num(result) - 1;
end

AddBarFiles(fragCountFiles, ...
            outputFile, ...
            sampleInterval, ...
            genomeChromSizes, ...
            1000000/numFragments);

