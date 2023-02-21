clear 
close all

addpath(genpath('/Users/jack/Documents/Github/1DVarS/test_cases/FRF2017/data'))

t = ncread('FRF-ocean_waves_8m-array_201611.nc','time');
t = t-t(1);
dt = (t(2) - t(1))/3600; % calculating time step in hours
date = 18;
hour = 16;
t_idx = 24*(date-1)+(hour-1);

wave_dir = ncread('FRF-ocean_waves_8m-array_201611.nc','waveDirectionBins');
dir_wed = ncread('FRF-ocean_waves_8m-array_201611.nc','directionalWaveEnergyDensity');
wave_f = ncread('FRF-ocean_waves_8m-array_201611.nc','waveFrequency');

figure()
s = pcolor(wave_dir,wave_f,dir_wed(:,:,t_idx)');
set(s,'EdgeColor','none');
title('Directional Wave Energy Density - 11/18/2016')
xlabel('wave angle (degrees)')
ylabel('wave frequency (Hz)')
