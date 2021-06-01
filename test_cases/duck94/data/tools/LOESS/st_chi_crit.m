function [crit_chi_sq] = st_chi_crit(m, prob)

% [crit_chi_sq] = st_chi_crit(m, prob)
%
% calculate crit_chi_sq such that 
%   prob[any_chi_sq <= crit_chi_sq] = prob
%
% input
%	m is degrees of freedom
%	prob, is probability distribution value
% output
%	crit_chi_sq is the chi-squared value giving 
%		specified probability


% proceed by bisection
% seed an itteration scheme with mean (= dof):
seed = m;
% and initial probability 
	F = gammainc(seed/2,m/2);

% and some limits for the itteration scheme
max_x = 6*m;   % this ought to be way out there
min_x = 0;

% cut the difference between the current seed and the maximum 
% seed in the correct direction this thing converges in a hurry, 
% so set a very small tollerance

while (abs(prob-F) > 0.00001)
	if ((F-prob)>0),
		% don't want to go higher
		max_x = seed;
		seed = (seed + min_x)/2;
		end
	if ((F-prob)<0),
		% don't want to go lower
		min_x = seed;
		seed = (seed + max_x)/2;
		end
	F = gammainc(seed/2,m/2);
end
crit_chi_sq = seed;
