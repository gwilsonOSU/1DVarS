function [zi, rmse, bi, ci, sk, Ni] = ...
   		loessInterp(x,z,sigma_z,xi,lx);
%
% 	function [zi, rmse, bi, ci, sk, Ni] = 
%			loessInterp(x,z,sigma_z,xi,lx);
%
% Filter/interpolate data z(x) to xi using a quadratic loess filter 
%	characterized by length scales lx; 
% 	(uses method of Schlax and Chelton, 1992)
%
% Input
%  x, the nxm observation location matrix, m is the number of dimensions
%  z, the nx1 observation array
%  sigma_z are the expected rms of observation errors on z
%     sigma_z could be constant (equal weight) or variable
%     If no information of this, input sigma_z = 1 (safe deffult choice)
%  xi, the pxm interpolation location matrix
%  lx are correlation scales, dx1
%      ok, I did synthetic tests with white spectrum, which suggest that:
%      -Quadratic filter: (6 dx) < lx < (Lpreserve/2)
%        where dx is sample spacing ...
%        and Lpreserve is wavelength that suffers less than 10% loss in variance
%
% Output
%    zi, are the estimated values at xi
%    rmse, are the estimated rms errors for the estimate
%    bi, are the estimated model coefficients
%    ci, are the estimated model coefficients rmse
%    sk, are the skill estimates of the weighted regression
%    Ni, are the number of effective dof at each location

%  Nathaniel Plant, 11/00

% measure input
[n,m] = size(x);
[ni, mi] = size(xi);
if (mi~=m)
   error('dimension of location matrices differ');
end

% fix up sigma_z, if constant
[ns,ms]=size(sigma_z);
if(max([ns,ms])==1)
   sigma_z = repmat(sigma_z,n,1);
else
   sigma_z = sigma_z(:);
end

% deal with nans
id = find(isfinite([z,x,sigma_z]*ones(m+2,1)));
n = length(id);
x=x(id,:);
z = z(id);
sigma_z = sigma_z(id);

% deal with entry of zero for obs. err.
max_sigma = max(sigma_z);
min_sigma = min(sigma_z);
if (min_sigma==0)
   % give the zero entries a weight equal to 1 order of magnitude greater 
   % than the smallest non-zero entry
   sigma_z(sigma_z==0) = sigma_z(sigma_z==0) + min(sigma_z(sigma_z>0))*1e-1;
   warning('Found sigma_z=0 and replaced with 0.1 min(sigma_z>0)')
end

% find out if any parameters are constant
xmean = mean(x);
xsdev = std(x);
idp = find(xsdev>0);
d=length(idp);

% center and normalize field area
x = x(:,idp)-ones(n,1)*xmean(idp);
x = x(:,idp)./(ones(n,1)*xsdev(idp));
xi = xi(:,idp)-ones(ni,1)*xmean(idp);
xi = xi(:,idp)./(ones(ni,1)*xsdev(idp));
lx = lx(:);
lx = lx(idp)./xsdev(idp)';

% remove global linear (for stability) trend
zmean = mean(z);
zsdev = std(z);
z = z - zmean;
XX = (x'*x)/n;
XZ = (z'*x)/n;
XX_inv = svdInvert(XX);
btrend = XX_inv*XZ'; 
ztrend = x*btrend;
zi_trend = zmean + xi*btrend;

% compute deviations from trend
zdev = z - ztrend;

% compute global variance
mso_dev = zdev'*zdev/n;
fprintf('standard deviation of obs. is %.3f\n', sqrt(mso_dev))

% pass each xi to model
zi = nan*xi(:,1);
rmse = zi;
sk = zi;
Ni = zi;
bi = nan*ones(ni,m);
ci = nan*ones(ni,m);
kx = 1./lx(:);
tic
for i = 1:ni
   % center on xi
   if (~mod(i,100))
	tnow = toc;
	temp_str = [num2str(i) ' of ' num2str(ni) '  ' num2str(tnow*(ni-i)/i) ...
			' seconds to go'];
	fprintf('\r%s     ', temp_str);
   end
   t = x-ones(n,1)*xi(i,:);
    % don't send all data
   id = find((abs(t)*kx)<d);

   if(length(id)>0)
      % interpolate one
      [zi(i), rmse(i), b, c, sk(i), Ni(i)] = ...
         loessFilt(t(id,:),zdev(id),sigma_z(id),lx);
      bi(i,idp) = (b'+ btrend')./xsdev(idp);
      ci(i,idp) = c'./xsdev(idp);
   end
end
zi = zi + zi_trend;
fprintf('\n');
%
% Copyright by Oregon State University, 2002
% Developed through collaborative effort of the Argus Users Group
% For official use by the Argus Users Group or other licensed activities.
%
% $Id: loessInterp.m,v 1.1 2004/08/20 00:51:06 stanley Exp $
%
% $Log: loessInterp.m,v $
% Revision 1.1  2004/08/20 00:51:06  stanley
% Initial revision
%
%
%key interp 
%comment  Quadratic loess interpolation master routine 
%
