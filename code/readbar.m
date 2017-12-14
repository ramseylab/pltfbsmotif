% reads .bar file to matlab
% [SEQNAME,data]=readbar(file)

% quick hack based on affy documentation
% Matti Nykter 20071105
% more hacks to make it faster 20080925

function [SEQNAME,data]=readbar(file)
fid=fopen(file,'r','b');

magic=char(fread(fid,8,'uchar')');
version=fread(fid,1,'float32');
NSEQ=fread(fid,1,'uint32');
NCOL=fread(fid,1,'uint32');
for i=1:NCOL
  types{i}=fread(fid,1,'uint32');
end
noparam=fread(fid,1,'uint32');
for i=1:noparam
  pnl=fread(fid,1,'uint32');
  param{i}.name=char(fread(fid,pnl,'uchar')');
  pnv=fread(fid,1,'uint32');
  param{i}.value=char(fread(fid,pnv,'uchar')');
end

for k=1:NSEQ
  SEQNAMELEN=fread(fid,1,'uint32');
  SEQNAME{k}=char(fread(fid,SEQNAMELEN,'uchar')');
  SEQGROUPNAMELEN=fread(fid,1,'uint32');
  SEQGROUPNAME=char(fread(fid,SEQGROUPNAMELEN,'uchar')');
  SEQVERLEN=fread(fid,1,'uint32');
  SEQVER=char(fread(fid,SEQVERLEN,'uchar')');
  noparvals=fread(fid,1,'uint32');
  for i=1:noparvals  
    pnl=fread(fid,1,'uint32');
    SEQparam{i}.name=char(fread(fid,pnl,'uchar')');
    pnv=fread(fid,1,'uint32');
    SEQparam{i}.value=char(fread(fid,pnv,'uchar')');
  end
  
  nodatapoints=fread(fid,1,'uint32');
  %datapoint=zeros(nodatapoints,NCOL);
  %for i=1:nodatapoints
  %  for j=1:NCOL
  %    if types{j}==2
  %      datapoint(i,j)=fread(fid,1,'int32');
  %    elseif types{j}==1
  %      datapoint(i,j)=fread(fid,1,'float32');
  %    else
  %      error('unknown datatype')
  %    end
  %  end
  %end
  if NCOL==2
    % to make it faster, read all data in large blocks (twice :)
    % read "first column"
    dp1=fread(fid,nodatapoints*NCOL,'int32');
    % rewind fid
    status=fseek(fid,-nodatapoints*NCOL*4,0);
    % read "second column"
    dp2=fread(fid,nodatapoints*NCOL,'float32');  
    % dump "wrong" values from both reads
    data{k}=[dp1(1:2:end),dp2(2:2:end)];
  else
     error('Error: fix the code!')
  end
end

fclose(fid);
