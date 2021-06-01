function out=smoothConditional(a,b,out_a_gt_b,out_a_lt_b,transition_range)
%
% out = smoothConditional(a,b,out_a_gt_b,out_a_lt_b,transition_range)
%
% A smooth differentiable function that approximates a simple if-else
% conditional statement.  The statement must be of the simple form
%
%   if(a<b)
%    out=out_a_gt_b;
%   else
%    out=out_a_lt_b;
%   end
%
% The output of smoothConditional() replaces the above if-else statement
% with a smoothly varying cosine-shaped function, ramping over an interval
% abs(a-b) < transition_range.
%
% Why do this?  Because the van der A model contains several non-smooth
% piecewise functions which causes discontinuous transport predictions and
% wreaks havoc with TL-AD code.
%

wgt=.5*(1+sin(pi*(a-b)/transition_range));
wgt(a-b>+transition_range/2)=1;
wgt(a-b<-transition_range/2)=0;
out = out_a_lt_b.*(1-wgt) + out_a_gt_b.*wgt;
