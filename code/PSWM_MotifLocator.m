function W = PSWM_MotifLocator(S,M,B,offset,W1,W2)

% W = PSWM_MotifLocator(S,M,B,offset,W1,W2) - Scan sequences using MotifLocator
%
% This function implements the traditional MotifLocator sequence scanner 
% using a position specific weight matrix (PSWM) M and a background model 
% B. All sequences in a matrix S (rows) are scanned using MotifLocator (and
% using an offset of the promoter start site). If precomputed results are
% specified in W1 for the MotifScanner and/or in W2 for the background,
% then they are used.
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
% B       - A (4^k)-by-4 matrix specifying a k:th order Markovian 
%           background model. Matrix consists of positive values such 
%           that each row sums up to one. M(i,j) defines the frequency of 
%           base pair j given that the previous k base pairs are i (as 
%           encoded into an integer) in the background B. An example: if
%           k=3, then i=1 corresponds to AAA, i=2 to AAC, i=3 to AAG, i=4
%           to AAT, i=5 to ACA, ..., i=64 to TTT.
% offset  - An optional parameter specifying the offset of the promoter 
%           start site. If offset>0, then base pairs S(:,1:offset) are 
%           ignored during matrix scanning (but can be used for a higher-
%           order (Markovian) background model).
% W1      - A matrix of size N-by-(L-l+1) consisting of the PSWM scores 
%           for all base pair positions and for all sequences, or empty
%           matrix.
% W2      - A matrix of size N-by-L consisting of the background model 
%           scores for all base pair positions and for all sequences, or
%           empty matrix.
%
% OUTPUT:
% W       - A matrix of size N-by-(L-l+1) consisting of the MotifLocator 
%           scores for all base pair positions and for all sequences.

% 12/10/2006 by Harri Lähdesmäki

[N,L] = size(S);
if nargin<4
  offset = round(log(size(B,1))/log(4)); % Background model order.
end
L = L - offset;
l = size(M,2);

if isempty(W1)
  W1 = PSWM_Scan(S,M,offset);
end

if isempty(W2)
  W2 = BM_Scan(S,B,offset);
end

W3 = zeros(N,L-l+1);
for i=1:L-l+1
  W3(:,i) = prod(W2(:,[i:i+l-1]),2);
end

W = W1./W3;
%W = W1./W2(:,1:L-l+1);
