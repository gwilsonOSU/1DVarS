function m_inv = svdInvert2(m);
% invert matrix m
% 
% m_inv = svdInvert(m);
%
% input
%    m is a NOT necessarily a square matrix
% output
%    m_inv is its inverse, in the SVD sense

%  Code from Nathaniel Plant, modified a tad 11/00 by Holman

SMALL_RATIO = 1e-30;		% define a small ratio
[Nt, Nx] = size(m);

% invoke svd
[u,d,v] = svd(m,0);

% Inversion is done by inverting the diagonal elements of d.  Elements that
% are too small (defined by ratio to d(1,1), since matlab orders eirenvalues)
% are neglected.

dinv = 1./diag(d);
not_small = find(dinv/dinv(1) < 1/SMALL_RATIO);
n = Nx-length(not_small);
dd = diag([dinv(not_small); zeros(n,1)]);
m_inv = v*dd'*u';		% now do the inverse thing

%if (n<Nx)
    %warning(sprintf('%d degenerate columns\n', n))
%end
%
% Copyright by Oregon State University, 2002
% Developed through collaborative effort of the Argus Users Group
% For official use by the Argus Users Group or other licensed activities.
%
% $Id: svdInvert.m,v 1.1 2004/08/20 00:51:07 stanley Exp $
%
% $Log: svdInvert.m,v $
% Revision 1.1  2004/08/20 00:51:07  stanley
% Initial revision
%
%
%key interp 
%comment  Matrix inversion by SVD 
%
