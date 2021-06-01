function ad_h=ad_calcTanbeta(ad_tanbeta,x)

nx=length(x);

ad_h=zeros(nx,1);

%3 tl_tanbeta(nx)=-(tl_h(nx)-tl_h(nx-1))/(x(nx)-x(nx-1));
ad_h(nx)  =ad_h(nx)  - 1./(x(nx)-x(nx-1)).*ad_tanbeta(nx);
ad_h(nx-1)=ad_h(nx-1)+ 1./(x(nx)-x(nx-1)).*ad_tanbeta(nx);
ad_tanbeta(nx)=0;

%2 tl_tanbeta(1)=-(tl_h(2)-tl_h(1))/(x(2)-x(1));
ad_h(2)=ad_h(2)- 1/(x(2)-x(1))*ad_tanbeta(1);
ad_h(1)=ad_h(1)+ 1/(x(2)-x(1))*ad_tanbeta(1);
ad_tanbeta(1)=0;

%1 tl_tanbeta(2:(nx-1))=-(tl_h(3:nx)-tl_h(1:(nx-2)))./(x(3:nx)-x(1:(nx-2)));
ad_h(3:nx)    =ad_h(3:nx)    - 1./(x(3:nx)-x(1:(nx-2))).*ad_tanbeta(2:(nx-1));
ad_h(1:(nx-2))=ad_h(1:(nx-2))+ 1./(x(3:nx)-x(1:(nx-2))).*ad_tanbeta(2:(nx-1));
ad_tanbeta(2:(nx-1))=0;
