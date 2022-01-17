function [A,dfddksd,dfdddelta,dfwdahat,dfwdksw] = ...
    vanderA_shields_Jacobian(d50,udelta,uhat,delta,mu,eta,lambda,ahat,...
                             ksd,ksw,fd,fw,...
                             branch_A1,branch_A4,branch_A5)

physicalConstants;

% eqn A4
if(branch_A4==1)
  % fw = .00251*exp(5.21*(ahat/ksw)^-.19);
  dfwdksw = .00251*exp(5.21*(ahat/ksw)^-.19) ...
            *5.21*(-.19)*(ahat/ksw)^(-.19-1)*(-ahat/ksw^2);
  dfwdahat=.00251*exp(5.21*(ahat/ksw)^-.19) ...
           *5.21*(-.19)*(ahat/ksw)^(-.19-1)/ksw;
else
  % fw = 0.3;
  dfwdksw=0;
  dfwdahat=0;
end

% eqn 20
% fd = 2*.4^2/log(30*delta/ksd)^2;
dfddksd=2*.4^2*(-2)/log(30*delta/ksd)^3 ...
        *(ksd/30/delta)*(-30*delta/ksd^2);
dfdddelta = 2*.4^2*(-2)/log(30*delta/ksd)^3 ...
    *ksd/30/delta*(30/ksd);

% eqn A1
if(branch_A1==1)
  % F=@(ksd,ksw)3*d90 + .4*eta^2/lambda;
  dFdksd=0;
  dFddelta=0;
  dFdksw=0;
  dFdd50=0;
  dFdd90=3;
  dFdmu=0;
  dFdudelta=0;
  dFduhat=0;
  dFdahat=0;
  if(lambda*eta==0)
    dFdeta=0;
    dFdlambda=0;
  else    
    dFdeta=2*.4*eta/lambda;
    dFdlambda=-.4*eta^2/lambda^2;
  end
else
  % F = @(ksd,ksw) d50*(mu-6) ...
  %     + 3*udelta^2/((s-1)*g)*fd(ksd) ...
  %     + 1.5*uhat^2/((s-1)*g)*fw(ksw,ahat) ...
  %     + .4*eta^2/lambda;
  dFdksd = 3*udelta^2/((s-1)*g)*dfddksd;
  dFddelta = 3*udelta^2/((s-1)*g)*dfdddelta;
  dFdksw = 1.5*uhat^2/((s-1)*g)*dfwdksw;
  dFdd50=mu-6;
  dFdd90=0;
  dFdmu=d50;
  dFdudelta=2*3*udelta/((s-1)*g)*fd;
  dFduhat = 2*1.5*uhat/((s-1)*g)*fw;
  dFdahat=1.5*uhat^2/((s-1)*g)*dfwdahat;
  if(lambda*eta==0)
    dFdeta=0;
    dFdlambda=0;
  else    
    dFdeta=2*.4*eta/lambda;
    dFdlambda=-.4*eta^2/lambda^2;
  end
end
% ksd = @(ksd,ksw) F(ksd,ksw);

% eqn A5
if(branch_A5==1)
  % G=@(ksd,ksw)d50 + .4*eta^2/lambda;
  dGdksd = 0;
  dGddelta = 0;
  dGdksw = 0;
  dGdahat = 0;
  dGdd50 = 1;
  dGdd90=0;
  dGdmu = 0;
  dGdudelta = 0;
  dGduhat = 0;
  if(lambda*eta==0)
    dGdeta=0;
    dGdlambda=0;
  else    
    dGdeta    = .4*2*eta/lambda;
    dGdlambda = -.4*eta^2/lambda^2;
  end
else
  % G=@(ksd,ksw) d50*(mu-6) ...
  %   + 3*udelta^2/((s-1)*g)*fd(ksd) ...
  %   + 1.5*uhat^2/((s-1)*g)*fw(ksw,ahat) ...
  %   + .4*eta^2/lambda;
  dGdksd = 3*udelta^2/((s-1)*g)*dfddksd;
  dGddelta = 3*udelta^2/((s-1)*g)*dfdddelta;
  dGdksw = 1.5*uhat^2/((s-1)*g)*dfwdksw;
  dGdd50=mu-6;
  dGdd90=0;
  dGdmu=d50;
  dGdudelta=2*3*udelta/((s-1)*g)*fd;
  dGduhat = 2*1.5*uhat/((s-1)*g)*fw;
  dGdahat=1.5*uhat^2/((s-1)*g)*dfwdahat;
  if(lambda*eta==0)
    dGdeta=0;
    dGdlambda=0;
  else    
    dGdeta=2*.4*eta/lambda;
    dGdlambda=-.4*eta^2/lambda^2;
  end
end
% ksw = @(ksd,ksw) G(ksd,ksw);

B = [dFdd50 dFdd90 dFdahat dFdmu dFduhat dFdudelta dFdeta dFdlambda dFddelta;
     dGdd50 dGdd90 dGdahat dGdmu dGduhat dGdudelta dGdeta dGdlambda dGddelta];
A = inv(eye(2) - [dFdksd dFdksw;
                  dGdksd dGdksw])*B;
