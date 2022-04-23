function tl_ubar=tl_dubarbierUmod(tl_ubar0,tl_kabs,tl_lambda,ubar0,kabs,lambda,x)

nx=length(x);

xb=lambda*2*pi./kabs;
tl_xb = -lambda*2*pi./kabs.^2.*tl_kabs ...
        + tl_lambda*2*pi./kabs;
for i=1:nx
  ind=find(x(i)-xb(i)<=x&x<=x(i));
  if(length(ind)<=1)
    tl_ur(i,1) = tl_ubar0(i,1);
    tl_ur(i,2) = tl_ubar0(i,2);
  else
    xx=xb(i)-x(i)+x(ind);
    xx=xx(:);
    for n=1:length(ind)
      tl_xx(n)=tl_xb(i);
    end
    tl_xx=tl_xx(:);
    term2=sum(xx);
    tl_term2=0;
    for n=1:length(ind)
      tl_term2 = tl_term2 + tl_xx(n);
    end
    for j=1:2
      term1=sum(xx.*ubar0(ind,j));
      tl_term1=0;
      for n=1:length(ind)
        tl_term1 = tl_term1 + tl_xx(n)*ubar0(ind(n),j) + xx(n)*tl_ubar0(ind(n),j);
      end
      tl_ur(i,j) = tl_term1/term2 - term1/term2^2*tl_term2;
    end
  end
end
tl_ubar = tl_ur;
