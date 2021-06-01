function tl_tanbeta=tl_calcTanbeta(tl_h,x)

nx=length(x);
tl_tanbeta(2:(nx-1))=-(tl_h(3:nx)-tl_h(1:(nx-2)))./(x(3:nx)-x(1:(nx-2)));
tl_tanbeta(1)=-(tl_h(2)-tl_h(1))/(x(2)-x(1));
tl_tanbeta(nx)=-(tl_h(nx)-tl_h(nx-1))/(x(nx)-x(nx-1));
