%
% assuming input mat-files exist in this directory, run the assim_1dh.m code
% to generate and save output mat-files
%
clear

load input_jul.mat
[posterior,diagnostics]=assim_1dh(prior,obs,1);
save assim_1dh_output_jul.mat

load input_oct.mat
[posterior,diagnostics]=assim_1dh(prior,obs,1);
save assim_1dh_output_oct.mat
