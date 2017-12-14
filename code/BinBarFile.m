function BinBarFile(barFileName, gcs, newBinSize, operator, outBarFileName)
% BinBarFile - Downsamples an Affymetrix BAR file using mean, median, etc.
%
%  BinBarFile(barFileName, gcs, newBinSize, operator,  outBarFileName)
%
% This function reads in numeric data track along the chromosome in
% the binary Affymetrix BAR file format, and subsamples the data
% along the chromosome in intervals of fixed size newBinSize.  The
% subsampling is done using a user-specified function, either min,
% max, median, or mean.  The resulting file is saved to a BAR file.
%
% Inputs:
%   barFileName - the name of the input BAR file
%   gcs - a struct array of chromosome names & sizes, made by the
%     function GetGenomeChromSizes.m
%   newBinSize - the sampling interval (in bp); a typical value 
%     would be 100.
%   operator - a string variable giving the name of the function
%     to apply to all the samples in an interval; must be one of
%     the following:  min, max, mean, median
%   outBarFileName - the name of the output BAR file to be written
%
% Uses: 
%   CellToWriteBar.m

% 
% =================================================
% BinBarFile.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================

[seqname, data]=readbar(barFileName);

C = length(seqname);

bucketCoords = cell(1,C);
buckets = cell(1,C);

for c=1:C
  chromName = seqname{c};
  chromInd = find(strcmp({gcs.chromName}, chromName));
  if isempty(chromInd)
    error('unable to find size for chromosome: %s', chromName);
  end
  chromSize = gcs(chromInd).chromSize;
  chromBucketCoords = 0:newBinSize:(chromSize-1);
  bucketCoords{c}=chromBucketCoords;
  B = length(chromBucketCoords);
  chromBuckets = zeros(1,B);
  dataCoords = data{c}(:,1);
  dataVals = data{c}(:,2);
  
  % compute the bucket indices of all the samples along the chromosome
  dataInds = min(B, max(1, 1 + floor((dataCoords+0.5*newBinSize)/newBinSize)));

  if strcmp(operator, 'mean')
    funcUse = @mean;
  else 
    if strcmp(operator, 'min') 
      funcUse = @min;
    else 
      if strcmp(operator, 'max')
        funcUse = @max;
      else
        if strcmp(operator, 'median')
          funcUse = @median;
        else
          error('unknown operator');
        end
      end
    end
  end
  
  [uniqueDataInds, junk, indPtrs]=unique(dataInds);
  chromBuckets(uniqueDataInds) = accumarray(indPtrs, dataVals, [length(uniqueDataInds) 1], funcUse);
  nzInds = find(chromBuckets);
  buckets{c} = [chromBucketCoords(nzInds)' chromBuckets(nzInds)'];
end

CellToWriteBar(outBarFileName, seqname, buckets);


