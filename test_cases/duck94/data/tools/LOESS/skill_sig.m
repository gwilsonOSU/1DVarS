function [skill_crit, pr] = skill_sig(skill,m,dof,alpha)
% [skill_crit, pr] = skill_sig(skill,m,dof,alpha)
%
% test if skill or correlation be significant at (1-alpha)% level
%
% input:
%	skill (or squared crosscorrelation) value
%   m is the number of model parameters to get skill=R2, pxy^2
%     2 is the correct number for linear model (slope+intercept are parameters)
%	dof, the effective number of degrees of freedom
%      get these from one of either long lag or other methods
% output:
%	skill_crit is the critical value above which we reject
%      hypothesus that true skill is zero
%	pr is probability at which skill (or P_xy^2) is significant

% chi_crit is critical one side chi_sq value for
% 95% confidence and two dof is 5.99 --> calculate it
% for m dof calculate it
chi_2_crit = st_chi_crit(m,1-alpha);

skill_crit = chi_2_crit./(dof + chi_2_crit);

% the thing that is chi_sq dist is:
X = dof*skill./(1-skill);

% probability that we would exceed this value 
pr = 1-st_prob_chisq(m, X);
