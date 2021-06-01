function t=Uwave_ruessink2012_tcrest(omega,r,phi,tguess)
%
% t=Uwave_ruessink2012_tcrest(omega,r,phi,tguess)
%
% Helper function for qtrans_vanderA.m: locate the time of maximum velocity
% in the crest (or trough), by setting du/dt = 0.  The initial guess
% 'tguess' should be close to a point where du/dt = 0, and this alone
% determines whether you are finding the time of the crest or trough.

dudt=@(t)cos(omega*t) ...
     - (sin(omega*t)+r*sin(phi)/(1+sqrt(1-r^2)))/(1-r*cos(omega*t+phi)) ...
     *r*sin(omega*t+phi);
t=abs(fzero(@(t)dudt(abs(t)),tguess));
