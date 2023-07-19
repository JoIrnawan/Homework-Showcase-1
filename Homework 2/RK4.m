function [output] = RK4(initialx,initialy,parameterh,finalx,function_handle)

% Given dy/dx = f(x,y)
% RK4 will numerically solve a first order ODE and output a set of x and y.

y(1) = initialy;   x(1) = initialx;
while max(x) < finalx
    k1 = function_handle(x(end),y(end));
    k2 = function_handle(x(end) + parameterh/2, y(end) + k1*parameterh/2);
    k3 = function_handle(x(end) + parameterh/2, y(end) + k2*parameterh/2);
    k4 = function_handle(x(end) + parameterh, y(end) + k3*parameterh);
    x(end + 1) = x(end) + parameterh;
    y(end + 1) = y(end) + (k1 + 2*k2 + 2*k3 + k4)*parameterh/6;
end
[output] = [x;y];
end