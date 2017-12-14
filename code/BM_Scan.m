function W = BM_Scan(S,B,offset)

% W = BM_Scan(S,B,offset)
%
% This function implements the traditional sequence scanning using a 
% higher-order Markovian background model. All sequences in a matrix S 
% (rows) are scanned with a Markovian model B of order k (using an offset 
% of the promoter start site).
%
% INPUT:
% S       - A matrix of size N-by-(L+offset) consisting of integers 
%           between 1 and 4. N is the number of sequences, L denotes 
%           the length of the sequences, and offset specifies the 
%           number of base pairs that are used to extend the sequences 
%           from the beginning/left. Thus, S(i,offset+1:offset+L) encodes 
%           the i:th sequence. Base pairs are coded using integers as 
%           follows: A = 1, C = 2, G = 3, and T = 4.
% B       - A (4^k)-by-4 matrix specifying a k:th order Markovian 
%           background model. Matrix consists of positive values such 
%           that each row sums up to one. M(i,j) defines the frequency of 
%           base pair j given that the previous k base pairs are i (as 
%           encoded into an integer) in the background B. An example: if
%           k=3, then i=1 corresponds to AAA, i=2 to AAC, i=3 to AAG, i=4
%           to AAT, i=5 to ACA, ..., i=64 to TTT.
% offset  - A parameter specifying the offset of the promoter start site. 
%           Parameter offset should satisfy offset>=k. If offset>k, then 
%           base pairs S(:,1:offset-k) are ignored during the scanning.%
%
% OUTPUT:
% W       - A matrix of size N-by-L consisting of the background model 
%           scores for all base pair positions and for all sequences.

% 10/10/2006 by Harri Lähdesmäki

[Bd1,Bd2] = size(B);
k = round(log(Bd1)/log(4)); % Background model order.

if nargin<3
  offset = k;
end

if offset<k
  error('offset needs to exceed the model order.');
end

[N,L] = size(S);
L = L - offset;

if k==0
  % In the case of 0:th ordel background model, use the standard 
  % PSWM scan.
  W = PSWM_Scan(S,B',offset);
  
else
  
  % Initialize the output.
  W = ones(N,L);
  
  offset_1 = offset - 1;
  offset_k = offset - k;
  b = 4.^[k-1:-1:0]';
  SS = S - 1; % Values between 0 and 3.
  
  for i=1:L
    ind = sub2ind([Bd1,Bd2],SS(:,offset_k+i:offset_1+i)*b+1,S(:,offset+i));
    W(:,i) = W(:,i).*B(ind);
  end
end
