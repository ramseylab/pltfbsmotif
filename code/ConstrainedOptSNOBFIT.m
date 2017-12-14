function xbest=ConstrainedOptSNOBFIT(fcn, u, v, x0, ncall, fglob, ftol, outputFileName)

prt = 0;
file = outputFileName;  % the intermediate data (after each call to SNOBFIT) 
                % are stored in test.mat
% fcn in {'bra','cam','gpr','hm3','hm6','ros','sh10','sh5','sh7',shu'}

n = length(u);  % dimension of the problem
% the following are meaningful default values
npoint = ceil(ncall/20);   % number of random start points to be generated
nreq = 4*n;     % no. of points to be generated in each call to SNOBFIT
x = rand(npoint,n);
T = (n+1)/3;
if T ~= ceil(T)
  error 'T should be an integer';
end

x(1:(npoint-1),:) = x(1:(npoint-1),:)*diag(v-u) + ones((npoint-1),1)*u'; 
x(1:(npoint-1),(T+1):2*T) = x(1:(npoint-1),1:T) + ...
                            rand(npoint-1,T).^0.1 .* ...
                            bsxfun(@minus, v((T+1):2*T)', x(1:(npoint-1),1:T));
x((npoint),:)=x0;

                % starting points in [u,v]
dx = (v-u)'*1.e-4; % resolution vector
p = 0.25;        % probability of generating a point of class 4
prt = 0;        % print level
                % prt = 0 prints ncall, xbest and fbest if xbest has
                %         changed
                % prt = 1 prints
	
'initial search of the parameter space'		
tic;
for j=1:npoint
  f(j,:) = [fcn(x(j,:)) sqrt(eps)];
  fprintf(1, '.');
% computation of the function values (if necessary, with additive
% noise)
end
[fbest fi] = min(f(:,1));
xbest = x(fi,:);
iterTime=toc/npoint;
ncall0 = npoint;   % function call counter
params = struct('bounds',{u,v},'nreq',nreq,'p',p); % input structure
% repeated calls to Snobfit
while ncall0 < ncall % repeat till ncall function values are reached
  timeToGo = (ncall - ncall0)*iterTime / 60;
  sprintf('func calls: %d; best val so far: %f; minutes to go: %f\n',ncall0,fbest,timeToGo)
  % (if the stopping criterion is not fulfilled first)
  if ncall0 == npoint  % initial call
    [request,xbest2,fbest2] = snobfit(file,x,f,params,dx);
%    ncall0,xbest,fbest
  else                 % continuation call
    [request,xbest2,fbest2] = snobfit(file,x,f,params);
  end
  if fbest2 < fbest
    fbest = fbest2;
    xbest = xbest2;
  end
  %  'xbest:'
%  for zz=1:n
%    disp(xbest(zz))
%  end
  if prt>0, request, end
  clear x
  clear f
  for j=1:size(request,1)
    x(j,:) = request(j,1:n);
    f(j,:) = [fcn(x(j,:)) sqrt(eps)];
    fprintf(1, '.');
  end 
  ncall0 = ncall0 + size(f,1); % update function call counter
  [fbestn,jbest] = min(f(:,1)); % best function value
  if fbestn < fbest
    fbest = fbestn;
    xbest = x(jbest,:);
%    'xbest:'
%    for zz=1:n
%      disp(xbest(zz))
%    end    
    ncall0,fbest % display current number of function values,
                       % best point and function value if fbest has
                       % changed
  end
  % check stopping criterion 
  % if fglob == 0, stop if fbest < 1.e-5
  % otherwise, stop if (fbest-fglob)/abs(fglob) < 1.e-2
  if fglob 
    if abs((fbest-fglob)/fglob) < ftol,break,end
  else
    if abs(fbest) < ftol,break,end
  end
end
%'xbest:'
for zz=1:n
%  disp(xbest(zz))
end
%ncall0,fbest  % show number of function values, best point and
	      % function value