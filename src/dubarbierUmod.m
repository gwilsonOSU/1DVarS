function ubar=dubarbierUmod(ubar0,kabs,lambda,x);

nx=length(x);

xb=lambda*2*pi./kabs;
for i=1:nx
  ind=find(x(i)-xb(i)<=x&x<=x(i));
  if(length(ind)<=1)
    ur(i,1)=ubar0(i,1);
    ur(i,2)=ubar0(i,2);
  else
    xx=xb(i)-x(i)+x(ind);
    xx=xx(:);
    term2=0;
    for n=1:length(ind)
      term2=term2+xx(n);
    end
    for j=1:2
      term1=0;
      for n=1:length(ind)
        term1 = term1 + xx(n)*ubar0(ind(n),j);
      end
      ur(i,j)=term1/term2;
    end
  end
end
ubar=ur;
