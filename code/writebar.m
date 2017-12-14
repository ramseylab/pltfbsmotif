% [chromo,chr_loc,redundancy]=textread('CoordStart','%s%d%d','headerlines',1);
% writebar('bartest.bar',properchrnames(chromo),chr_loc,log2(data))
%

% quick hack based on affy documentation
% Matti Nykter 20071105

% speed up
% Theo Knijnenburg Sep 30 2009

function writebar(file,chromo,chr_loc,data)
%[chromo,chr_loc,redundancy]=textread('CoordStart','%s%d%d','headerlines',1);
%file='bartest.bar'
fid=fopen(file,'W','b');

magic='barr1234';
fwrite(fid,magic,'uchar');
version=2;
fwrite(fid,version,'float32');

SEQNAMES=unique(chromo);
NSEQ=length(SEQNAMES);
fwrite(fid,NSEQ,'uint32');
NCOL=2;
fwrite(fid,NCOL,'uint32');
types=[2,1];
for i=1:NCOL
  fwrite(fid,types(i),'uint32');
end

noparam=0;
fwrite(fid,noparam,'uint32');
param={'file_type','scale'};
paramv={'intensity','log2'};
for i=1:noparam
  pnl=length(param{i});
  fwrite(fid,pnl,'uint32');
  fwrite(fid,param{i},'uchar');
  pnv=length(paramv{i});
  fwrite(fid,pnv,'uint32');
  fwrite(fid,paramv{i},'uchar');
end

for k=1:NSEQ
  SEQNAMELEN=length(SEQNAMES{k});
  fwrite(fid,SEQNAMELEN,'uint32');
  SEQNAME=[SEQNAMES{k}];
  fwrite(fid,SEQNAME,'uchar');
  SEQGROUPNAMELEN=0;
  fwrite(fid,SEQGROUPNAMELEN,'uint32');
  SEQGROUPNAME='';
  fwrite(fid,SEQGROUPNAME,'uchar');
  SEQVERLEN=0;
  fwrite(fid,SEQVERLEN,'uint32');
  SEQVER='';
  fwrite(fid,SEQVER,'uchar');
  
  noparvals=0;
  fwrite(fid,noparvals,'uint32');
  SEQparam={};
  SEQparamv={};
  for i=1:noparvals  
    pnl=length(SEQparam{i});
    fwrite(fid,pnl,'uint32');
    fwrite(fid,SEQparam{i},'uchar');
    pnv=length(SEQparamv{i});
    fwrite(fid,pnv,'uint32');
    fwrite(fid,SEQparamv{i},'uchar');
  end

  idx=find(strcmp(SEQNAMES{k},chromo)==1);
  nodatapoints=length(idx);
  fwrite(fid,nodatapoints,'uint32');

  %before speed up
%   for i=1:nodatapoints
%     for j=1:NCOL
%       if types(j)==2
%         fwrite(fid,chr_loc(idx(i)),'int32');
%       elseif types(j)==1
%         fwrite(fid,data(idx(i)),'float32');
%       else
%         error('unknown type')
%       end
%     end
%   end

%   %after speed up 1 --- seemd to crash every now and then
  fwrite(fid,chr_loc(idx(1)),'int32');
  fwrite(fid,chr_loc(idx(2:nodatapoints)),'int32',4);
  fseek(fid,-nodatapoints*8+8,0);
  fwrite(fid,data(idx(1)),'float32');
  fwrite(fid,data(idx(2:nodatapoints)),'float32',4);
  
  %after speed up 2 (In batches of N datapoints)
%   N = 1e5;
%   vec = [1:N:nodatapoints nodatapoints+1];
%   for i=1:length(vec)-1
%       beginpos  = vec(i);
%       endpos    = vec(i+1)-1;
%       lengthpos = endpos-beginpos+1;
%       fwrite(fid,chr_loc(idx(beginpos)),'int32');
%       fwrite(fid,chr_loc(idx(beginpos+1:endpos)),'int32',4);
% %       position = ftell(fid)
%       fseek(fid,-lengthpos*8+8,0);
% %       position = ftell(fid)
%       fwrite(fid,data(idx(beginpos)),'float32');
%       fwrite(fid,data(idx(beginpos+1:endpos)),'float32',4);
%   end
  
  
end

fclose(fid);
