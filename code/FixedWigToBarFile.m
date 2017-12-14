function FixedWigToBarFile(wigFile, barFile)

[status, result]=system(['cat ' wigFile ' | grep -v fixedStep | wc -l | cut -f1']);
N = str2num(result) - 1;
fid=fopen(wigFile, 'r');
notAtEnd = 1;
lastChrom = '';
chroms = cell(N,1);
coords = zeros(N,1);
values = zeros(N,1);
lineCtr = 0;
while 1
  tline = fgets(fid);
  if tline == -1
    break;
  end
  if strncmp(tline, 'track type', 10)
    continue;
  end
  tokens=regexp(tline, '^fixedStep chrom=(\S+) start=(\d+) step=(\d+)', ...
                'tokens', ...
                'once');
  if ~isempty(tokens)
    chrom = tokens{1};
    % WIG files use coordinates numbered from 1, so we need to
    % subtract 1 from the chromosomal coordinate
    coord = str2num(tokens{2}) - 1;
    stepSize = str2num(tokens{3});
    stepCtr = 0;
  else
    value = str2num(tline);
    stepCtr = stepCtr + 1;
    lineCtr = lineCtr + 1;
    chroms{lineCtr} = chrom;
    realCoord = coord + (stepSize * (stepCtr - 1)) - 1;
    coords(lineCtr) = realCoord;
    values(lineCtr) = value; 
  end
end

fclose(fid);

writebar(barFile, chroms, coords, values);

