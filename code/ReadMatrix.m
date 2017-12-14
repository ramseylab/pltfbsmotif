function m=ReadMatrix(matFile)
% Reads a motif matrix file in clover format:
%
% >V$MYOD_01
% 1 2 2 0
% 2 1 2 0
% 3 0 1 1
% 0 5 0 0
% 5 0 0 0
% 0 0 4 1
% 0 1 4 0
% 0 0 0 5
% 0 0 5 0
% 0 1 2 2
% 0 2 0 3
% 1 0 3 1

[as, cs, gs, ts]=textread(matFile, '%f %f %f %f', 'headerlines', 1);
N = length(as);
m = zeros(4, N);
m(1,:)=as;
m(2,:)=cs;
m(3,:)=gs;
m(4,:)=ts;
for n=1:N
  m(:,n) = m(:,n) / sum(m(:,n));
end
