function [sequencePartInds]=SplitSequenceByNs(seqNums, margin)
L = length(seqNums);
nanInds = find(isnan(seqNums));
partCtr = 0;
sequencePartInds = {};
lastStart = 0;
N = length(nanInds);
if N > 0
  if nanInds(1) > margin
    sequencePartInds{1} = 1:(nanInds(1)-1);
    partCtr = partCtr + 1;
  end
  for n=2:N
    if nanInds(n) > nanInds(n-1) + margin 
      % we have skipped a scannable sequence
      partCtr = partCtr + 1;
      sequencePartInds{partCtr} = (nanInds(n-1)+1):(nanInds(n)-1);
    end
  end
  if L > nanInds(N) + margin
    sequencePartInds{partCtr+1} = (nanInds(N)+1):L;
  end
else
  sequencePartInds{1} = 1:L;
end


