function seqScores=ScanSequenceForMatrix(sequence, matrix, background)

if ischar(sequence)
  seqNums = basepairs2num(sequence);
else
  seqNums = sequence;
end

L = length(sequence);
M = size(matrix, 2);
seqPartIndsList=SplitSequenceByNs(seqNums, M); 
P = length(seqPartIndsList);
seqScores = nan(1,L);
k = round(log(size(background,1))/log(4));
for p=1:P
  seqPartInds = seqPartIndsList{p};
  seqPart = seqNums(seqPartInds);
  seqPartScoresFwd = PSWM_MotifLocator(seqPart, matrix, background, k, [], []);
  seqPartScoresRev = fliplr(PSWM_MotifLocator(oppositestrand(seqPart), matrix, background, k, [], []));

  seqScores(seqPartInds((k+1):(length(seqPartInds)-M+1))) = max(seqPartScoresFwd, ...
                                                            seqPartScoresRev);
end


