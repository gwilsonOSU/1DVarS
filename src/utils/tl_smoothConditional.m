function tl_out=tl_smoothConditional(tl_a,tl_b,tl_out_a_gt_b,tl_out_a_lt_b,a,b,out_a_gt_b,out_a_lt_b,transition_range)
%
% TL code for smoothConditional.m
%

wgt=.5*(1+sin(pi*(a-b)/transition_range));
tl_wgt = .5*cos(pi*(a-b)/transition_range).*pi*( ...
    + (tl_a-tl_b)/transition_range ...
    - (a-b)/transition_range.^2.*tl_transition_range );

wgt(a-b>+transition_range/2)=1;
tl_wgt(a-b>+transition_range/2)=0;

wgt(a-b<-transition_range/2)=0;
tl_wgt(a-b<-transition_range/2)=0;

tl_out = tl_out_a_lt_b.*(1-wgt) ...
         - out_a_lt_b.*tl_wgt ...
         + tl_out_a_gt_b.*wgt ...
         + out_a_gt_b.*tl_wgt;
