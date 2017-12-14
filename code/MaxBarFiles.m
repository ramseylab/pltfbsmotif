function MaxBarFiles(inputBarFiles, outputBarFile)

I = length(inputBarFiles);

[chromNames, data]=readbar(inputBarFiles{1});
C = length(chromNames);
chromSizes = zeros(1,C);
dataAllFiles = cell(1,C);

outData = cell(1,C);
for c=1:C
  chromInData = data{c};
  chromSize = size(chromInData,1);
  chromSizes(c) = chromSize;
  dataAllFiles{c} = zeros(chromSize, I);
  chromOutData = zeros(chromSize, 2);
  chromOutData(:,1) = chromInData(:,1);
  outData{c} = chromOutData;
end

barData = cell(1,I);

for i=1:I
  inputBarFiles{i}
  if i > 1
    [chromNames, data]=readbar(inputBarFiles{i});
    if length(chromNames) ~= C
      error 'invalid bar file';
    end
  end
  for c=1:C
    chromSize = chromSizes(c);
    chromDataAllFiles = dataAllFiles{c};
    chromInputData = data{c};
    chromDataAllFiles(:,i) = chromInputData(:,2);
    dataAllFiles{c} = chromDataAllFiles;
  end
end

for c=1:C
  chromDataAllFiles = dataAllFiles{c};
  chromDataMax = max(chromDataAllFiles, [], 2);
  chromOutData = outData{c};
  chromOutData(:,2) = chromDataMax;
  outData{c} = chromOutData;
end

CellToWriteBar(outputBarFile, chromNames, outData);
