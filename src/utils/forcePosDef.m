function C = forcePosDef(C)
%
% C = forcePosDef(C)
%
% Force a matrix to be positive definite by eliminating negative
% eigenvalues.  If you are using this it might be because you have a bug in
% the code that's screwing up your matrices.  But, it can be a useful hammer.

C = .5*(C + C');  % symmetric
[V,D]=eig(C);
D=diag(D);
D(D<0)=0;  % enforce +ve def
C=V*diag(D)*V';