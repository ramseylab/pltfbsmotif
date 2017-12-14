function numFragments = SumExtendedFragmentFiles(fragCountFiles, fragFilesDir)

F = length(fragCountFiles);

numFragments = 0;
for f=1:F
  fileName = fragCountFiles{f};
  [pathstr, name]=fileparts(fileName);
  bedFileName=[fragFilesDir '/' name '.bed'];
  [status, result]=system(['cat ' bedFileName ' | wc -l | cut -f1']);
  numFragmentsFile = str2num(result)
  if isempty(numFragmentsFile) || 0==numFragmentsFile
    error('unable to find extended fragment file: %s', bedFileName);
  end
  numFragments = numFragments + numFragmentsFile - 1;
end
