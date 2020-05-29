function [x,y] = circle_xy(x0,y0,r,n)
% returns x and y coordinates of circle centered at [x0,y0] with radius r
%
% optional forth argument specifies number of points around the circle
% (default is 101)

if nargin < 4
    n = 101;
end

th = linspace(0,2*pi,n);
x = r * cos(th) + x0;
y = r * sin(th) + y0;