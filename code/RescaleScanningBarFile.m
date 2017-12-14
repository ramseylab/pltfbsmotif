function RescaleScanningBarFile(inputBarFile, outputBarFile, motifLength)

RescaleBarFile(inputBarFile, outputBarFile, @(x) DoRescaling(x, motifLength));

end

function outputData=DoRescaling(inputData, motifLength)
outputData = ones(size(inputData));
indPos = find(inputData > 1);
outputData(setxor(1:length(inputData), indPos)) = 0;
outputData(indPos) = log10(inputData(indPos));
outputData(indPos) = outputData(indPos) - min(outputData(indPos));
outputData(indPos) = outputData(indPos) / median(outputData(indPos));
end
