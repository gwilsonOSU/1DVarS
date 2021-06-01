function t = findBarAndDeep(t, bar, d)
%   trans = findBarAndDeep(trans, bar, deep)
%
%  find intersection points of a transect with the bar position
%  data and with the deep data

% find the equations of line segments forming the offshore bar.
x = [bar.xs]; y = [bar.ys];
m = diff(y)./diff(x);
m(isinf(m)) = 10^3;         % 89.94 degrees, close enough.
b = y(1:end-1) - m.*x(1:end-1);

% find intersection of transect with bar segments.  Chose intersection
% that is closest to each transect line
yInt = (b - t.b*m/t.m) ./ (1-m/t.m);
yInt(t.m==0) = t.xy0(2);         % fix segments with zero slope
pick = find((yInt>=(y(1:end-1))) & (yInt<=(y(2:end))));
if isempty(pick)
    if yInt>y(end)  % beyond the end 
        pick = length(y)-1;
    else
        pick = 1;
    end
end
pick = pick(1);         % to avoid multiple options
t.bar.x = (yInt(pick)-b(pick))/m(pick);
t.bar.y = yInt(pick);
t.bar.d = sqrt((t.bar.x-t.xy0(1))^2+(t.bar.y-t.xy0(2))^2); % along trans distance to bar

% finally, find intersections of transects and deep segments.  
% Again, chose intersection that is closest to each transect line

% find the equations of line segments forming the offshore deep contour.
x = [d.x]; y = [d.y];h = [d.h];beta = [d.beta];
m = diff(y)./diff(x);
m(isinf(m)) = 10^2;         % 89.4 degrees, close enough.
b = y(1:end-1) - m.*x(1:end-1);

yInt = (b - t.b*m/t.m) ./ (1-m/t.m);
yInt(t.m==0) = t.xy0(2);
pick = find((yInt>(y(1:end-1))) & (yInt<(y(2:end))));
if isempty(pick)
    if yInt>y(end)
        pick = length(y)-1;
    else
        pick = 1;
    end
end
pick = pick(1);         % some chance of having two choices
t.deep.x = (yInt(pick)-b(pick))/m(pick);
t.deep.y = yInt(pick);
t.deep.d = sqrt((t.deep.x-t.xy0(1))^2+ ...
            (t.deep.y-t.xy0(2))^2); % along trans distance to deep
t.deep.h = h(pick) + (h(pick+1)-h(pick))/(y(pick+1)-y(pick))* ...
    (t.deep.y-y(pick));
t.deep.beta = beta(pick) + (beta(pick+1)-beta(pick))/(y(pick+1)-y(pick))* ...
    (t.deep.y-y(pick));
    