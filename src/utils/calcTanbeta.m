function tanbeta=calcTanbeta(x,h)

nx=length(x);
tanbeta(2:(nx-1))=-(h(3:nx)-h(1:(nx-2)))./(x(3:nx)-x(1:(nx-2)));
tanbeta(1)=-(h(2)-h(1))/(x(2)-x(1));
tanbeta(nx)=-(h(nx)-h(nx-1))/(x(nx)-x(nx-1));
