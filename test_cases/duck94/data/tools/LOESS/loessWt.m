function w = loessWt(x,lx)
%
% [w] = loessWt(x,lx)
%
% Input
%   x are nxm inputs, weights will be centered on x=0
%   lx is mx1 length scales
%
% Output
%   w are weights 0<=w<=1

[n,m] = size(x);
lx = lx(:);

% normalize by correlation scale
x = x./(ones(n,1)*(lx'));

% weights for least-sq. minimization
x = x.^2;

% sum across all inputs, if m>1
if(m>1)
   x = x*ones(m,1);
end

% the loess weighting function (zero if x>1)
w = ((1-(x.^3)).^3).*(x<1);  % used by greenslade
% w = ((1-(x.^1.5)).^3).*(x<1); % used by schlax and chelton
%
% Copyright by Oregon State University, 2002
% Developed through collaborative effort of the Argus Users Group
% For official use by the Argus Users Group or other licensed activities.
%
% $Id: loessWt.m,v 1.1 2004/08/20 00:51:07 stanley Exp $
%
% $Log: loessWt.m,v $
% Revision 1.1  2004/08/20 00:51:07  stanley
% Initial revision
%
%
%key interp 
%comment  Calculates weights for weighted regression 
%
