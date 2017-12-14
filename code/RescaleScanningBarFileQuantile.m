function RescaleScanningBarFileQuantile(inputBarFile, outputBarFile)

RescaleBarFile(inputBarFile, outputBarFile, @(x) DoRescaling(x));

end

function outputData=DoRescaling(inputData)
  outputData = tiedrank(inputData)/length(inputData);
end
