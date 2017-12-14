function CombineTracksIntoMatrix(trackFiles, outputMatrixFile, normalizeToOne)

T = length(trackFiles);

[chromNames, data]=readbar(trackFiles{1});

C = length(chromNames);
F = 0;
for c=1:C
  F = F + size(data{c},1);
end

chromInds = zeros(F,1);
chromCoords = zeros(F,1);

trackData = zeros(F,T);

for t=1:T
  trackFiles{t}
  if t > 1
    [seqNames, data]=readbar(trackFiles{t});
  end
  featCtr = 1;
  if normalizeToOne(t)
    maxMax = -inf;
    for c=1:C
      chromData = data{c};
      maxMax = max(maxMax, max(chromData(:,2)));
    end
  end
  for c=1:C
    chromData = data{c};
    if nargin > 2
      if normalizeToOne(t)
        chromData(:,2) = chromData(:,2)/maxMax;
      end
    end
    chromFeats = size(chromData,1);
    trackData(featCtr:(featCtr + chromFeats - 1),t) = chromData(:,2);
    if t==1
      chromCoords(featCtr:(featCtr + chromFeats - 1)) = chromData(:,1);
      chromName = chromNames{c};
      chromInd = find(strcmp(chromNames, chromName));
      if isempty(chromInd)
        error 'could not find chromosome';
      end
      chromInds(featCtr:(featCtr + chromFeats - 1)) = chromInd;
    end
    featCtr = featCtr + chromFeats;
  end
end

save(outputMatrixFile, ...
     'trackData', ...
     'chromNames', ...
     'chromInds', ...
     'chromCoords', ...
     'trackFiles');
