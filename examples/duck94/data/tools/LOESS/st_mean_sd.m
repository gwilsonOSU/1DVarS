function [u, s, n] = st_mean_sd(data, bad_flag)
% [u, s, n] = st_mean_sd(data, bad_flag)
%
% calculate mean and sd of data
%
% input
%   data, input data series (1d)
%   bad_flag, flags bad or missing data
% output:
%	u is data mean
%	s is sd
%	n is number of good points

% used nan as default flag
if (nargin==1)
   bad_flag=nan;
end

% check for misaligned vector
N = size(data);
if (N(1)==1 | N(2)==1)
   data = data(:);
   N(2) = 1;
elseif (N(1)==0)
   u=bad_flag;
   s=bad_flag;
   n=0;
   return;
end

% initialize output
u=bad_flag*ones(1,N(2));
s=bad_flag*ones(1,N(2));
n=bad_flag*ones(1,N(2));
for i=1:N(2)
   id = find(data(:,i)~=bad_flag & isfinite(data(:,i)));
   n(i)=length(id);
   if (n(i)>0)
      u(i) = mean(data(id,i));
      s(i) = sqrt(mean((data(id,i)-u(i)).^2));
   end
end