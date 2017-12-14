function fileSum=SumBarFile(barFileName)
%SumBarFile - sum up all the data samples in an Affy BAR file
%
%  fileSum=SumBarFile(barFileName)
%
% Computes the number of the samples in an Affymetrix BAR file 
% whose values exceed a user-specified threshold value. Returns
% this number, as well as the total number of samples in the
% BAR file.
% 
% Inputs: 
%   barFile - the full name of the BAR file to be analyzed
%   threshold - the sample value threshold
%
% Outputs:
%   total - the total number of samples in the BAR file
%   above - the number of samples in the BAR file that are
%           above the threshold value
%
% Subfunctions:            none
% MAT-files required:      none
% See also:                none
%

% =================================================
% BarFileCountThresh.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================

if nargin < 2
  threshold = 0;
end

[seqNames, data]=readbar(barFileName);
C = length(seqNames);
fileSum = 0;
for c=1:C
  vals = data{c}(:,2);
  fileSum = fileSum + sum(vals);
end



