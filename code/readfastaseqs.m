function [S,L] = readfastaseqs(file)

% [S,L] = readfastaseqs(file)
%
% This function reads in fasta sequences from a fasta file.
%
% INPUT:
% file  - Name of a fasta file.
%
% OUTPUT:
% S     - Cell vector of sequences.
% L     - Cell vector of sequence labels.


fid = fopen(file);
S = {};
L = {};

str = '';
i = 0;
while 1 % Read the file line by line.
  
  l = fgetl(fid); % The current line.
  
  if isempty(l) % Skip empty lines.
    continue;
    
  elseif isnumeric(l) % Break if end of file (l==-1)
    if i>0
      S{end+1} = str; % Store the previous sequence.
    end
    break;
    
  elseif l(1)=='>' % Start of a new sequence
    i = i + 1;
    if i>1
      S{end+1} = str; % Store the previous sequence.
    end
    if length(l)==1
      error('No sequnce name.');
    end
    L{end+1} = l(2:end);
    str = '';
    
  else
    str = [str,l];
  end
  
end

fclose(fid);
