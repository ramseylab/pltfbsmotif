ypred = [ 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 ];

bindingSites = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];

ypreddiff = diff(ypred);

predInd = find(ypred);

transInds = find(ypreddiff);

% a false positive occurs when "ypred" has a 1, and "bindingSites"
% has a zero, in the same slot; thus we want to look at the values
% of "bindingSites" at all slots where "ypred" has a 1:

bindingAtPredSites = bindingSites(predInd);

% but we don't want to just count up the zero elements of this
% vector, because that would overcount the number of false
% positives.  If the code makes two adjacent positive predictions,
% it shouldn't be penalized.

