function W = PSWM_Scan(S,M,offset)

% W = PSWM_Scan(S,M,offset) - Scan a set of sequences
%
% This function implements the traditional sequence scanning using a 
% position specific weight matrix (PSWM). All sequences in a matrix S 
% (rows) are scanned with a PSWM matrix M (using an offset of the promoter 
% start site).
%
% INPUT:
% S       - A matrix of size N-by-(L+offset) consisting of integers 
%           between 1 and 4. N is the number of sequences, L denotes 
%           the length of the sequences, and offset specifies the 
%           number of base pairs that are used to extend the sequences 
%           from the beginning/left. Thus, S(i,offset+1:offset+L) encodes 
%           the i:th sequence. Base pairs are coded using integers as 
%           follows: A = 1, C = 2, G = 3, and T = 4.
% M       - A 4-by-l matrix specifying a position specific weight matrix 
%           (PSWM). Matrix consists of positive values such that each
%           column sums up to one. M(i,j) defines the frequency of base
%           pair i at position j in the site/motif M.
% offset  - An optional parameter specifying the offset of the promoter 
%           start site. If offset>0, then base pairs S(:,1:offset) are 
%           ignored during the scanning (but can be used for a higher-
%           order (Markovian) background model.
%
% OUTPUT:
% W       - A matrix of size N-by-(L-l+1) consisting of the PSWM scores 
%           for all base pair positions and for all sequences.

% 10/10/2006 by Harri Lähdesmäki

if nargin<3
  offset = 0;
end

[N,L] = size(S);
L = L - offset;
l = size(M,2);

% Initialize the output.
W = ones(N,L-l+1);

offset_1 = offset - 1;

for i=1:(L-l+1)
  for j=1:l
    W(:,i) = W(:,i).*M(S(:,offset_1+i+j),j);
  end
end
