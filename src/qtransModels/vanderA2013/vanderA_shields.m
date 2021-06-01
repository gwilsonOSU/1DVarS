function [theta_av,ksd,ksw,fd,fw,...
          branch_A1,branch_A4,branch_A5] = ...
    vanderA_shields(d50,d90,udelta,uhat,delta,mu,eta,lambda,ahat)
%
% [theta_av,ksd,ksw,fd,fw,...
%  branch_A1,branch_A4,branch_A5] = ...
%    vanderA_shields(d50,d90,udelta,uhat,delta,mu,eta,lambda,ahat)
%
% Solves the nonlinear system of 6 eqns defined by vanderA's (A1)-(A5) and
% (20).  Is constructed in such a way as to facilitate TL-AD coding, by
% setting up the problem as a nonlinear system, and keeping track of which
% branches of max/min and piecewise functions are used in the solution.
%
% NOTE: it is assumed for convenience 'udelta' is actually the magnitude of
% the mean current, not a vector.  This just makes it easier than writing
% 'udeltaabs' all throughout the code
%

physicalConstants;

% hack: if lambda=0 or eta=0, then just set lambda=1 and eta=0.  The point is just to
% ignore ripples if they don't occur.  Will need to be careful about this in
% TL coding.
if(eta*lambda==0)
  lambda=1;
  eta=0;
end

% eqn A3
theta_av = @(fd,fw) (.5*fd*udelta^2 + .25*fw*uhat^2)/((s-1)*g*d50);

% eqn 20.  Note, see helper function below for fw, eqn A4
fd = @(ksd) 2*(.4/log(30*delta/ksd))^2;

% eqn A1
ksd = @(ksd,ksw) max( 3*d90, d50*(mu+6*(theta_av(fd(ksd),fw_fn(ksw,ahat))-1)) ) ...
      + .4*eta^2/lambda;

% eqn A5
ksw = @(ksd,ksw) max(   d50, d50*(mu+6*(theta_av(fd(ksd),fw_fn(ksw,ahat))-1)) ) ...
      + .4*eta^2/lambda;

% solve nonlinear system of equations for ksd, ksw
opt=optimset('Display','off');
kswguess = fzero(@(kk)ksw(abs(kk),abs(kk))-abs(kk),d50);  % first guess, assume ksd==ksw
kvecguess=[max(3*d90,kswguess); kswguess];
kvec = fsolve(@(kvec)[ksd(real(kvec(1)),real(kvec(2)))-real(kvec(1));
                      ksw(real(kvec(1)),real(kvec(2)))-real(kvec(2))],...
              kvecguess,opt);

% set all the output values, overriding the inline functions defined above
ksd=real(kvec(1));
ksw=real(kvec(2));
fd=fd(ksd);
fw=fw_fn(ksw,ahat);
theta_av=theta_av(fd,fw);

% take note of which branches were taken in the various max/min and
% piecewise functions used above.  This info is needed when defining the TL
% version of the nonlinear system.
if(3*d90 > d50*(mu+6*(theta_av-1)))
  branch_A1=1;
else
  branch_A1=2;
end
if(ahat/ksw>1.587)
  branch_A4=1;
else
  branch_A4=2;
end
if(d50 > d50*(mu+6*(theta_av-1)))
  branch_A5=1;
else
  branch_A5=2;
end

end  % end of main function

% eqn A4, written as helper function because of branch
function fw_out = fw_fn(ksw,ahat)
if(ahat/ksw>1.587)
  fw_out=.00251*exp(5.21*(ahat/ksw)^-.19);
else
  fw_out=0.3;
end
end  % end of helper function fw()
