function [output] = LinFit(x,y)

% Given the x and the y of a sample data values, LinFit evaluates the
% linear regression of the sample data.
% LinFit outputs the line of best fit with x values given.

X = [ones(length(x),1) x];
b = X\y;
y_fit = X*b;
output = [x,y_fit];
end