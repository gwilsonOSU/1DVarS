function [tl_theta_av,tl_ksd,tl_ksw,tl_fd,tl_fw] = ...
    tl_vanderA_shields(tl_d50,tl_d90,tl_udelta,tl_uhat,tl_delta,...
                       tl_mu,tl_eta,tl_lambda,tl_ahat,...
                       d50,d90,udelta,uhat,delta,mu,eta,lambda,ahat,...
                       ksd,ksw,fd,fw,...
                       branch_A1,branch_A4,branch_A5)

physicalConstants;

% Jacobian matrix
[A,dfddksd,dfdddelta,dfwdahat,dfwdksw] = ...
    vanderA_shields_Jacobian(d50,udelta,uhat,delta,mu,eta,lambda,ahat,...
                             ksd,ksw,fd,fw,...
                             branch_A1,branch_A4,branch_A5);

% solve nonlinear system of equations for ksd, ksw
% [ksd; ksw] = [F(ksd,ksw); G(ksd,ksw)];
% NL: uses fsolve()
% TL: use Jacobian matrix.  See notes in vanderA_shields_deriv/
tl_rhs = [tl_d50 tl_d90 tl_ahat tl_mu tl_uhat tl_udelta tl_eta tl_lambda tl_delta]';
tl_kvec = A*tl_rhs;  % use this line for AD calculation
tl_ksd=tl_kvec(1);
tl_ksw=tl_kvec(2);

% back-calculate other TL vars, given tl_ksd and tl_ksw
tl_fd = dfddksd*tl_ksd + dfdddelta*tl_delta;
tl_fw = dfwdahat*tl_ahat + dfwdksw*tl_ksw;
% theta_av = (.5*fd*udelta^2 + .25*fw*uhat^2)/((s-1)*g*d50);
tl_theta_av = -(.5*fd*udelta^2 + .25*fw*uhat^2)/((s-1)*g*d50^2)*tl_d50 ...
    + 1/((s-1)*g*d50)*( ...
        .5*tl_fd*udelta^2 ...
        + 2*.5*fd*udelta*tl_udelta ...
        + .25*tl_fw*uhat^2 ...
        + 2*.25*fw*uhat*tl_uhat );
