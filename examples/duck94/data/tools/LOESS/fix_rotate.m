function [xr, yr] = fix_rotate(x,y,theta,xo,yo);

% rotate some data
%
% [xr, yr] = fix_rotate(x,y,theta,xo,yo);
%
% Input
%   x,y are the un-rotated data 
%   theta is angle to rotate
%		(decimal degrees!)
%   xo,yo is the rotation origin
% 
% Output
%    xr, yr are rotated coords
%

% Example: FRF rotation problem (Kent's vs OSU (ngp))
% theta_kent = -18.14182902;
% theta_osu = ddmmss_dd(-18.0734) = -18.1261;
% dtheta = theta_kent - theta_osu =   -0.0157
% xo = 28.075;
% yo = 521.707;

% translate to origin
xr = x-xo;
yr = y-yo;

% get current azimuth
az = atan2(yr,xr); % radians

% get radius (meters)
r = sqrt((xr.^2) + (yr.^2));

% rotate
xr = r.*cos(az+theta);
yr = r.*sin(az+theta);

% translate away from origin
xr = xr+xo;
yr = yr+yo;

