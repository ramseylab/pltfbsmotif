function n = basepairs2num(s)

% n = basepairs2num(s)
%
% This function converts base pairs (4-mers) into integers 1,...,4.
% Use mapping: a = 1, c = 2, g = 3, and t = 4.
%
% INPUT:
% s   - A sequence/string of base pairs (a, c, g and t)
%
% OUTPUT:
% n   - Integer version of the sequence.
%
% Originally written by Harri Lahdesmaki, 2006
% Modified by Stephen Ramsey to accomodate N's, 2009

n = zeros(1,length(s));
s = lower(s);

for i=1:length(s)
  if s(i)=='a'
    n(i) = 1;
  elseif s(i)=='c'
    n(i) = 2;
  elseif s(i)=='g'
    n(i) = 3;
  elseif s(i)=='t'
    n(i) = 4;
  elseif s(i)=='n'
    n(i) = nan;
  else
    %warning('Unexpected base pair symbol.');
    %n = [];
    %return;
    n(i) = -1;
  end
end

if 0
  ind = find(n>0);
  if length(ind)<length(s)
    %warning('Unexpected base pair symbol ignored.');
    n = n(ind);
  end
end
