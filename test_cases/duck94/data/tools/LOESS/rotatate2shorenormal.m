function [x,y,theta_xy,xi,yi,XI,YI] = rotatate2shorenormal(x,y,z,e,xi,yi)
% [xr,yr,theta_xy,xir,yir,XIr,YIr] = rotatate2shorenormal(x,y,z,e,xi,yi)
%
% rotate data to a shore-normal coordinate where x is cross-shore and y is alongshore
%
% Input
%   x,y,z,e, the data (Nx1) each
%          e is the rms error on each obs
%   xi,yi, the target (Px1 or PxQ) each x,y coords (i.e., grid)
%
% Output
%   xr,yr, the rotated data coords
%   theta_xy, the rotation angle
%   xir,yir, the rotated target (grid) coords) 
%      e.g., [xr, yr] = fix_rotate(x,y,theta_xy,xrot0,yrot0);
%            where xrot0,yrot0 is center of rotation
%   XIr,YIr, a rotate coordinate regular grid 

% figure out the grid
[p,q] = size(xi);
if(p>1 & q>1)
    dxi = diff(xi(1,1:2));
    dyi = diff(yi(1:2,1));
    xi = xi(:);
    yi = yi(:);
end

% get center of rotation
xrot0 = mean(x);
yrot0 = mean(y);

% convert error to a weight and regress
w = e/std(z);
[beta_xy, beta_xy_ci, beta_xy_skill] = regr_xzw([ones(length(x),1),x,y], z, w);

% get rotation to downhill
theta_xy = atan2(beta_xy(3),-beta_xy(2));

% rotate data
[x, y] = fix_rotate(x,y,theta_xy,xrot0,yrot0); % data
% rotate grid
[xi, yi] = fix_rotate(xi,yi,theta_xy,xrot0,yrot0); % original grid

% for fast interpolation, we may want a grid in the local coordinate
XI = [floor(min(xi))]:dxi:[max(xi)]; % local regular grid in rotated frame
YI = [floor(min(yi))]:dyi:[max(yi)];
[XI, YI] = meshgrid(XI, YI); % NOTE, this means we always do 2-d interpolation
