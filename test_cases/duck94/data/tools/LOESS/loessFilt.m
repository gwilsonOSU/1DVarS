function [zo, rmse, b, c, sk, n] = loessFilt(X,z,sigma_z,lx);
%
%   function [zo, rmse, b, c, sk, n] = loessFilt(X,z,sigma_z,lx);
%
% return quadradic loess filter estimate for a single grid point at x=0
%
% Input
%   X are independent variables, usually [xs, ys], centered at prediction
%	location x=0.  i.e. [x0, y0] must be subtracted from the [xs, ys]
%	before passing the data to this routine. 
%       X is size nxd, where d is dimension (say 2 for x,y data or 3 for x,y,t.
%       Do NOT include constant in any columns of X
%   z are observations, nx1
% 	IMPORTANT, z are assumed to be somehow detrended.  
%   	Large mean values over the width, lx, confuse the error analysis
%   sigma_z are the expected observation errors on z
%      if errors unknown, then sigma_z=1 should be used
%   lx are correlation scales, dx1
%      ok, I did synthetic tests with white spectrum:
%      -Quadratic filter: (6 dx) < lx < (Lpreserve/2)
%      where dx is sample spacing 
%      and Lpreserve suffers less than 10% loss in variance
%
% Output
%    zo is the estimate at x=0
%    rmse is the estimated error
%    b are coefficients (excluding zo)
%    c are coefficients variances (excluding zo)
%    sk are model skill of regression
%    n are the effective dof for each regression

%  From Nathaniel Plant, 11/00

porder = 2;		% assume second order

% initialize output
[n,d] = size(X);
zo = nan;
rmse = nan;
sk=nan;
b = nan*ones(d,1);
c = nan*ones(d,1);

% get weights
w = loessWt(X,lx(:));

% should exit here if weights are zero
if(max(w)==0)
   warning('weights are all zero')
   return
end

% weight by variance
% don't worry about magnitude of weights, that is managed in the regression code
w = w./sigma_z(:);
w = w/max(w);

% build the input
X=[ones(n,1),X];
m=1+d;

% get quadratic terms
if(porder==2)
   for i = 2:(d+1)
      % get the quadratic self-terms
      m = m + 1;
      X(:,m) = (X(:,i).^2);
      for j = i:d
         % get the cross-terms
         if (j~=i)
            m = m + 1;
            X(:,m) = X(:,i).*X(:,j);
         end
      end
   end
end

% effective dof
nu = sum(w);

% check to see if we are overdetermined
if (nu<m)
   % fprintf('nu/m = %.1f / %d, not enough data\n', nu, m)
   return
end

% remove mean and send to regression, with weights
[b,c,sk, n, msz] = regrXZW(X,z,w);

% put the primary output in its place
zo = b(1);
rmse = c(1);

% get the 1st order terms
b = [b(2:(d+1))];
c = [c(2:(d+1))];

% now, if mse is negative
if (imag(rmse)~=0)
   rmse = nan;
   zo = nan;
   b = nan*b;
   c = b;
   return
end
%
% Copyright by Oregon State University, 2002
% Developed through collaborative effort of the Argus Users Group
% For official use by the Argus Users Group or other licensed activities.
%
% $Id: loessFilt.m,v 1.1 2004/08/20 00:51:06 stanley Exp $
%
% $Log: loessFilt.m,v $
% Revision 1.1  2004/08/20 00:51:06  stanley
% Initial revision
%
%
%key interp 
%comment  Quadratic loess interpolator for a point 
%
