function CompareMatrixToBinding(scanningBar, bindingBar, gcs)

[seqsScanning, dataScanning]=readbar(scanningBar);
[seqsBinding, dataBinding]=readbar(bindingBar);

C = length(gcs);

scanningValuesBinding = [];
scanningValuesNoBinding = [];


for c=1:C
  chromName = gcs(c).chromName;
  barInd = find(strcmp(seqsScanning, chromName));
  if ~isempty(barInd) 
    chromScanningData = dataScanning{barInd};
    chromBindingData = dataBinding{barInd};
    if size(chromScanningData,1) ~= size(chromBindingData,1)
      size(chromScanningData)
      size(chromBindingData)
      error 'inconsistent bar file sizes';
    end
    X = size(chromBindingData,1);
    indBinding = find(chromBindingData(:,2));
    scanningValuesBinding = [scanningValuesBinding; ...
                    chromScanningData(indBinding,2)];
    scanningValuesNoBinding = [scanningValuesNoBinding;
                    chromScanningData(setdiff(1:X, indBinding),2)];
  end
end

scanningValuesBinding = scanningValuesBinding(scanningValuesBinding > 0);
scanningValuesNoBinding = scanningValuesNoBinding(scanningValuesNoBinding > 0);

[h,p]=kstest2(scanningValuesBinding, scanningValuesNoBinding);
fprintf('P-value is: %e\n', p);
figure;
subplot(1,2,1);
hist(scanningValuesBinding);
title('With Binding');
subplot(1,2,2);
hist(scanningValuesNoBinding);
title('No Binding');

medBinding = median(scanningValuesBinding);
medNoBinding = median(scanningValuesNoBinding);
minVal = min([scanningValuesBinding; scanningValuesNoBinding]);

fprintf('median with binding: %f  no binding: %f  min value: %f\n', medBinding, medNoBinding, minVal);






