function nr = oppositestrand(n)

% nr = oppositestrand(n)
%
% This function converts a forward/reverse sequence strand into the 
% corresponding reverse complementary strand (or the other way around).
% Code assumes mapping: a = 1, c = 2, g = 3, and t = 4. Note that the
% returned sequence is also flipped around so that it can be directly used
% in scanning (with original matrices).
%
% INPUT:
% n   - A sequence
%
% OUTPUT:
% nr  - Reverse version of the sequence.

nr = zeros(1,length(n));
nr(n==1) = 4;
nr(n==2) = 3;
nr(n==3) = 2;
nr(n==4) = 1;
nr(isnan(n)) = nan;

nr = nr(end:-1:1);
