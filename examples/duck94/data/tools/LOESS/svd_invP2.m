function m_inv = svd_invP(m, p);
% invert matrix m
% retain eigen functions that account for p% of the variance
% 
% m_inv = svd_invP(m, p);
%
% input
%    m is a NOT necessarily a square matrix
%    p is the percent variance to retain
%
% output
%    m_inv is its inverse, in the SVD sense

% the retained variance (does not need to be 90%)
P = p/100;

[Nt, Nx] = size(m);

% invoke svd
[u,d,v] = svd(full(m));

% total variance
D = diag(d);
S = sum(D);

% dd is used in the inverse calculation
% if sum from 1 to i explains P*100% of variance, truncate
dd=zeros(Nt,Nx);
n=0;
for i = 1:Nx,
	if (sum(D(1:i))/S <= P)
		dd(i,i) = 1/d(i,i);
	else,
		n=n+1;
	end
end

if (n>0)
    %fprintf('%d un-used columns\n', n);
end

% now do the inverse thing
m_inv = v*dd'*u';
