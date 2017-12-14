function ConvertBarToWig(barFileName, wigFileName)
% ConvertBarToWig - Reads an Affymetrix BAR file and saves it in Wiggle format
%
%  ConvertBarToWig(barFileName, wigFileName)
%
% This function converts a numeric data track on one or more
% chromosomes, from the binary Affymetrix BAR file format to a text
% file in the UCSC Wiggle variableStep format.
%
% Inputs:
%  barFileName - the name of the Affymetrix BAR file to read
%  wigFileName - the name of the UCSC Wiggle file to write
%
% Uses: 
%  readbar.m (written by Matti Nykter)

% =================================================
% ConvertBarToWig.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================

% read the BAR file data
[seqname, data] = readbar(barFileName);

% get the number of chromosomes
C = length(seqname);

fid=fopen(wigFileName, 'w+');
fprintf(fid, 'track type=wiggle_0\n');
% go through each chromosome
for c=1:C 
  chromName=seqname{c};
  if strcmp(chromName(1:3),'chr')
    chromName=chromName(4:length(chromName));
  end
  fprintf(fid, 'variableStep chrom=%s\n', seqname{c});
  
  N = size(data{c},1);
  for n=1:N
    if 0==mod(n,100000)
      sprintf('For chromosome %s, have written %d out of %d data points', seqname{c}, n, N)
    end
    % WIG file coordinates are numbered starting at 1, so we need
    % to add 1 to the BAR file coordinate
    fprintf(fid, '%d\t%f\n', data{c}(n,1) + 1, data{c}(n,2));
  end
end
fclose(fid);


