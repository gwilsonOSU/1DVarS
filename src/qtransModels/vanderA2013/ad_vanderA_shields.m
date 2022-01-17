function [ad_d50,ad_d90,ad_udelta,ad_uhat,ad_delta,...
          ad_mu,ad_eta,ad_lambda,ad_ahat] = ...
    ad_vanderA_shields(ad_theta_av,ad_ksd,ad_ksw,ad_fd,ad_fw,...
                       d50,udelta,uhat,delta,mu,eta,lambda,...
                       ahat,ksd,ksw,fd,fw,...
                       branch_A1,branch_A4,branch_A5)

% Jacobian matrix
[A,dfddksd,dfdddelta,dfwdahat,dfwdksw] = ...
    vanderA_shields_Jacobian(d50,udelta,uhat,delta,mu,eta,lambda,ahat,...
                             ksd,ksw,fd,fw,...
                             branch_A1,branch_A4,branch_A5);

physicalConstants;

% init AD variables
ad_d50=0;
ad_d90=0;
ad_udelta=0;
ad_uhat=0;
ad_delta=0;
ad_mu=0;
ad_eta=0;
ad_lambda=0;
ad_ahat=0;

% NOTE: TL code supplies the matrix A s.t. tl_kvec = A*tl_rhs, and some
% other calculated derivatives supplied here as inputs

% %b2 back-calculate other TL vars, given tl_ksd and tl_ksw
%3 tl_theta_av = -(.5*fd*udelta^2 + .25*fw*uhat^2)/((s-1)*g*d50^2)*tl_d50 ...
%     + 1/((s-1)*g*d50)*( ...
%         .5*udelta^2*tl_fd ...
%         + 2*.5*fd*udelta*tl_udelta ...
%         + .25*uhat^2*tl_fw ...
%         + 2*.25*fw*uhat*tl_uhat );
ad_d50=ad_d50-(.5*fd*udelta^2 + .25*fw*uhat^2)/((s-1)*g*d50^2)*ad_theta_av;
ad_fd    =ad_fd     + 1/((s-1)*g*d50)*.5*udelta^2   *ad_theta_av;
ad_udelta=ad_udelta + 1/((s-1)*g*d50)*2*.5*fd*udelta*ad_theta_av;
ad_fw    =ad_fw     + 1/((s-1)*g*d50)*.25*uhat^2    *ad_theta_av;
ad_uhat  =ad_uhat   + 1/((s-1)*g*d50)*2*.25*fw*uhat *ad_theta_av;
ad_theta_av=0;
%2 tl_fw = dfwdahat*tl_ahat + dfwdksw*tl_ksw;
ad_ahat=ad_ahat+ dfwdahat*ad_fw;
ad_ksw =ad_ksw + dfwdksw *ad_fw;
ad_fw=0;
%1 tl_fd = dfddksd*tl_ksd + dfdddelta*tl_delta;
ad_ksd  =ad_ksd  +dfddksd  *ad_fd;
ad_delta=ad_delta+dfdddelta*ad_fd;
ad_fd=0;

% %b1 system of equations for ksd, ksw.
% [tl_ksd; tl_ksw] = A*[tl_d50 tl_d90 tl_ahat tl_mu tl_uhat tl_udelta tl_eta tl_lambda]';
% Note A is an 2x9 matrix.  After augmenting this TL system of eqns
% with the state variables ksd and ksw, the "full" TL model is a 11x11 matrix.
% The AD code corresponds to transposing this.
Ap=[[eye(9); A] zeros(11,2)];  % full 11x11 matrix for TL model
ad_rhs=[ad_d50 ad_d90 ad_ahat ad_mu ad_uhat ad_udelta ad_eta ad_lambda ad_delta ad_ksd ad_ksw]';
ad_rhs = Ap'*ad_rhs;
ad_d50   =ad_rhs(1);
ad_d90   =ad_rhs(2);
ad_ahat  =ad_rhs(3);
ad_mu    =ad_rhs(4);
ad_uhat  =ad_rhs(5);
ad_udelta=ad_rhs(6);
ad_eta   =ad_rhs(7);
ad_lambda=ad_rhs(8);
ad_delta =ad_rhs(9);
ad_ksd=0;
ad_ksw=0;
