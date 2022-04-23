function [ad_ubar0,ad_kabs,ad_lambda]=ad_dubarbierUmod(ad_ubar,ubar0,kabs,lambda,x)

nx=length(x);

% init AD
ad_ubar0=zeros(nx,2);
ad_kabs=zeros(nx,1);
ad_lambda=0;
ad_xb=zeros(nx,1);  % init
ad_ur=zeros(nx,2);   % init

xb=lambda*2*pi./kabs;
%3 tl_ubar = tl_ur;
ad_ur = ad_ur + ad_ubar;
% ad_dumpout.ubar=ad_ubar;
ad_ubar=0;
for i=1:nx
  ind=find(x(i)-xb(i)<=x&x<=x(i));
  if(length(ind)<=1)
    %2a2 tl_ur(i,2) = tl_ubar0(i,2);
    ad_ubar0(i,2)=ad_ubar0(i,2)+ad_ur(i,2);
    ad_ur(i,2)=0;
    %2a1 tl_ur(i,1) = tl_ubar0(i,1);
    ad_ubar0(i,1)=ad_ubar0(i,1)+ad_ur(i,1);
    ad_ur(i,1)=0;
  else
    xx=xb(i)-x(i)+x(ind);
    xx=xx(:);
    term2=sum(xx);
    ad_xx=zeros(length(ind),1);  % init
    ad_term2=0;  % init
    for j=1:2
      term1=sum(xx.*ubar0(ind,j));
      ad_term1=0;  % init
      %2b7 tl_ur(i,j) = tl_term1/term2 - term1/term2^2*tl_term2;
      ad_term1=ad_term1+ 1/term2      *ad_ur(i,j);
      ad_term2=ad_term2- term1/term2^2*ad_ur(i,j);
      ad_ur(i,j)=0;
      for n=1:length(ind)
        %2b6 tl_term1 = tl_term1 + tl_xx(n)*ubar0(ind(n),j) + xx(n)*tl_ubar0(ind(n),j);
        ad_xx(n)          =ad_xx(n)          + ubar0(ind(n),j)*ad_term1;
        ad_ubar0(ind(n),j)=ad_ubar0(ind(n),j)+ xx(n)          *ad_term1;
      end
      %2b5 tl_term1=0;
      ad_term1=0;
    end
    for n=1:length(ind)
      %2b4 tl_term2 = tl_term2 + tl_xx(n);
      ad_xx(n)=ad_xx(n)+ad_term2;
    end
    %2b3 tl_term2=0;
    ad_term2=0;
    %2b2 tl_xx=tl_xx(:);
    for n=1:length(ind)
      %2b1 tl_xx(n)=tl_xb(i);
      ad_xb(i)=ad_xb(i)+ad_xx(n);
      ad_xx(n)=0;
    end
  end
end
% tl_xb = -lambda*2*pi./kabs.^2.*tl_kabs ...
%         + tl_lambda*2*pi./kabs;
ad_kabs=ad_kabs- lambda*2*pi./kabs.^2.*ad_xb;
ad_lambda = ad_lambda + sum(2*pi./kabs.*ad_xb);
ad_xb=0;
