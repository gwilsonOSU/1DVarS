function [ad_a,ad_b,ad_out_a_gt_b,ad_out_a_lt_b]=ad_smoothConditional(ad_out,a,b,out_a_gt_b,out_a_lt_b,transition_range)
%
% AD code for smoothConditional.m
%

% init AD
nx=length(a.*b);
ad_a=zeros(nx,1);
ad_b=zeros(nx,1);
ad_out_a_gt_b=zeros(nx,1);
ad_out_a_lt_b=zeros(nx,1);

% NL bkgd
wgt=.5*(1+sin(pi*(a-b)/transition_range));
wgt(a-b>+transition_range/2)=1;
wgt(a-b<-transition_range/2)=0;

%4 tl_out = tl_out_a_lt_b.*(1-wgt) ...
%          - out_a_lt_b.*tl_wgt ...
%          + tl_out_a_gt_b.*wgt ...
%          + out_a_gt_b.*tl_wgt;
ad_out_a_lt_b=ad_out_a_lt_b+ (1-wgt)                  .*ad_out;
ad_out_a_gt_b=ad_out_a_gt_b+ wgt                      .*ad_out;
ad_wgt       =ad_wgt       + (out_a_gt_b - out_a_lt_b).*ad_out;
ad_out=0;

%3 tl_wgt(a-b<-transition_range/2)=0;
ad_wgt(a-b<-transition_range/2)=0;

%2 tl_wgt(a-b>+transition_range/2)=0;
ad_wgt(a-b>+transition_range/2)=0;

%1 tl_wgt = .5*cos(pi*(a-b)/transition_range).*pi*( ...
%     + (tl_a-tl_b)/transition_range ...
%     - (a-b)/transition_range.^2.*tl_transition_range );
coef=.5*cos(pi*(a-b)/transition_range).*pi*;
ad_a               =ad_a               + coef/transition_range          .*ad_wgt;
ad_b               =ad_b               - coef/transition_range          .*ad_wgt;
ad_transition_range=ad_transition_range- coef.*(a-b)/transition_range.^2.*ad_wgt;
ad_wgt=0;
