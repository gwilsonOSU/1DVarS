function [g,alpha,beta0,nu,rho,hmin,gammaType,betaType]=hydroParams();
%
% common params for NL-TL-AD model
%

g=9.8126;
alpha=1;
nu=.1;
rho=1030;
hmin=0.5;

% switch for breaker model.  Use '2001' Battjes and Stive (1985) (as used by
% Ruessink et al. 2001), or '2003' for Ruessink et al. (2003).
gammaType=2003;

% switch for roller model.  Use 'none' to turn off the roller, 'const' to
% set constant beta=beta0 for all gridpoints, or 'rafati21' to compute
% x-dependent beta based on Rafati et al., (2021) eqn 10, which they
% recommended for duck94 simulations
betaType='const';
beta0=0.1;
