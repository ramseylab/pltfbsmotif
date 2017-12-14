function RescaleBarFile(inputBarFile, outputBarFile, scalingFunc)

[seqs, data]=readbar(inputBarFile);

C = length(seqs);
for c=1:C
  chrData = data{c};
  chrData(:,2) = feval(scalingFunc, chrData(:,2));
  data{c} = chrData;
end

CellToWriteBar(outputBarFile, seqs, data);
