function BarToTSV(barFileName, tsvFileName, includeBarHeights)

[seqNames, data]=readbar(barFileName);

C = length(seqNames);

if nargin < 3
  includeBarHeights = 0;
end

fid=fopen(tsvFileName, 'w+');
for c=1:C
  chrData = data{c};
  F = size(chrData,1);
  for f=1:F
    if includeBarHeights
      fprintf(fid, '%s\t%d\t%f\n', seqNames{c}, chrData(f,1), chrData(f,2));
    else
      fprintf(fid, '%s\t%d\n', seqNames{c}, chrData(f,1));
    end
  end
end
fclose(fid);