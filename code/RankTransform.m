% RankTransformRows
%
% Stephen Ramsey, Institute for Systems Biology
% Mar. 2006
%
% Performs a rank transformation on the matrix "M" over the
% dimension "dim".  Each element in "rank" contains the rank
% of the corresponding element of "M" in the column (or row)
% of "M", based on whether "dim" is set to 1 or 2. 
%
function ranks=RankTransform(v)

[uniqueVals, m, n] = unique(v);
U = length(uniqueVals);
ranks = n/U;
