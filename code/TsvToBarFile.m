function TsvToBarFile(tsvFile, barFile)

[chrs, coords] = textread(tsvFile, '%s %d');

S = length(chrs);
uniqueChrs = unique(chrs);
C = length(uniqueChrs);
data = cell(1,C);
for c=1:C
  chrName = uniqueChrs{c};
  chrInds = find(strcmp(chrs, chrName));
  S = length(chrInds);
  chrData = ones(S,2);
  chrData(:,1) = sort(coords(chrInds));
  data{c} = chrData;
end

CellToWriteBar(barFile, uniqueChrs, data);

