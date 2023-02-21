function [F] = st_prob_chisq(m, input_chi_sq)
% [F] = st_prob_chisq(m, input_chi_sq)
%
% calculate probability that chi_sq <= input_chi_sq
% input
%    m is degrees of freedom
%    input_chi_sq is the chi_sq value that is being tested
% output
%    F is probability, fraction


% simple gamma function integral is evaluated by matlab:
F = gammainc(input_chi_sq/2,m/2);
