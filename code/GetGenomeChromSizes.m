function genomeChroms = GetGenomeChromSizes(genomeChromSizesFile)
%GetGenomeChromSizesFile - read genome chromosome sizes file
%
%  genomeChroms = GetGenomeChromSizes(genomeChromSizesFile)
%
% Returns a data structure containing the sizes of all 
% chromosomes in a genome.  The chromosome sizes are
% obtained from a tab-delimited file that should have
% been previously generated using MakeGenomeChromSizesFile.m.
%
% Inputs:
%   genomeChromSizesFile - the tab-delimited file 
%     into which the chromosome size data are saved
%
% Outputs:
%   genomeChroms - a struct array (length equal to the 
%     number of chromosomes in the genome), containing the
%     chromosome size data.  The fields are:
%        chromName - chromosome name
%        chromSize - chromosome size, in bp
%
% See also:    MakeGenomeChromSizes.m

% =================================================
% GetenomeChromSizesFile.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================
[genomeChromNames, genomeChromSizes] = textread(genomeChromSizesFile, ...
                                                '%s %d', ...
                                                'delimiter', '\t');

genomeChroms = struct('chromName', genomeChromNames, ...
                      'chromSize', num2cell(genomeChromSizes));
