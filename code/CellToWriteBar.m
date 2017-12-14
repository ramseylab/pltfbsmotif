function CellToWriteBar(barFileName, chromNames, data)
%CellToWriteBar - save a single genomic data track to a BAR file
%
%   CellToWriteBar(barFileName, chromNames, data)
%
% This is a convenience function that allows providing a 
% genomic data track in a cell array format, with each component
% of the cell array being the vector of data for a single
% chromosome, and to have this data written to a BAR file.
% The data are then reformatted to the format used by writebar.m.
% 
% Inputs:
%   barFileName - the name of the BAR file to be written
%   chromNames - a cell array (length C) of chromosome names
%   data - a cell array (length C), each entry of which contains
%       a two-column matrix of data.  The first column contains
%       the coordinate locations of the data track ("bars"), and
%       the second column contains the data values at those 
%       locations.
%
% Outputs: 
%
% Example:  
%
% Other m-files required:  writebar.m
% Subfunctions:            none
% MAT-files required:      none
% See also:                none
% 

% =================================================
% CellToWriteBar.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================
C = length(chromNames);
numPoints = 0;
for c=1:C
  coords = data{c}(:,1);
  numPoints = numPoints + length(coords);
end
coordsVec = zeros(numPoints,1);
dataVec = zeros(numPoints,1);
chromNamesInd = zeros(numPoints,1);
pointCtr = 1;
for c=1:C
  coords = data{c}(:,1);
  numPointsChr = length(coords);
  coordsVec(pointCtr:(pointCtr + numPointsChr - 1)) = coords;
  dataVec(pointCtr:(pointCtr + numPointsChr - 1)) = data{c}(:,2:end);
  chromNamesInd(pointCtr:(pointCtr + numPointsChr - 1)) = repmat(c, numPointsChr, 1);
  pointCtr = pointCtr + numPointsChr;
end

writebar(barFileName, ...
         chromNames(chromNamesInd), ...
         coordsVec, ...
         dataVec);

