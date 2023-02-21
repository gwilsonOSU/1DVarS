function [r,ar] = loess_kernelND(d,p,DX)

%fprintf('is it ok for d>1?\n')
% Input
%   d, the dimension of the kernel
%   p, the power of the loess (1=linear, 2=quadratic)
%
% Output
%   r, the radial distance
%   ar, the weights, sum to one, are scaled by number of inputs

if (nargin < 3 | isempty(DX)) % (~exist('DX'))
	DX = 0.05;
end

x = [-(1+DX):DX:(1+DX)]';
N = length(x);

% build ND input

switch d
case 1
    X = x;
case 2
    [X,Y] = ndgrid(x,x);
    X = [X(:),Y(:)];
case 3
    [X,Y,T] = ndgrid(x,x,x);
    X = [X(:),Y(:),T(:)];
end

% get weights

W = loess_wt(X);

% get basis for radially symmetric output
r = [x,zeros(N,d-1)];
w = loess_wt(r);

% build quadratic input
m = d+1;
n = N^d;

if(p==1)
   % build linear input
   X=[ones(n,1),X];
   r = [ones(N,1), r];
elseif(p==2)
   % get quadratic terms
   % there will be (m^2 + m)/2 entries (including the constant)

   q = 0.5*m*(m+1);
   X=[ones(n,1),X, zeros(n,q-m)];
   r=[ones(N,1),r, zeros(N,q-m)];

   for i = 2:(d+1)
      % get the quadratic self-terms
      for j = 2:i
         m  = m+1;
         X(:,m) = (X(:,i).*X(:,j));
         r(:,m) = (r(:,i).*r(:,j));
      end
   end
end

% weight all components
XW = [X.*repmat(W,1,m)];
r = [r.*repmat(w,1,m)];

% get the scaling matrix-- solves linear regression
XX = (XW'*XW)/(N^d);
XX_inv = inv(XX);

% we can get radially symmetric result now
% and this is the thing that is mult against the weighted data
ar = [r*XX_inv(:,1)].*w/N; 
r = x; 
return

% or continue to see the beast in N-dimensional space
% and this is the thing that is mult against the weighted data
bX = [XW*XX_inv(:,1)]; 

% so put the weight in this thing, rather than data, to get the indicator
% function
a = (bX(:,1).*W); 

% and reshape
switch d
case 1
    A = a;
otherwise
    A = reshape(a,N*ones(1,d));
end

% here is the 3-d surface in 1-d
ar = A(:,(N+1)/2, (N+1)/2);
