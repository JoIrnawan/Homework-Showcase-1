function [output] = RK4_2(initialx,initialy,initialDy,parameterh,finalx,D2y_handle,Dy_handle)

% Given D2y = f(x,y,Dy) and Dy = g(x,y,Dy) where D = d/dx
% 2RK4 will solve a second order ODE numerically using Runge-Kutta 4 and
% output a set of y, Dy, and x.

y(1) = initialy;   x(1) = initialx;   Dy(1) = initialDy;
while x(end) < finalx
    m1 = Dy_handle(x(end),y(end),Dy(end));
    k1 = D2y_handle(x(end),y(end),Dy(end));
    
    m2 = Dy_handle(x(end) + parameterh/2,...
        y(end) + parameterh*m1/2, Dy(end) + parameterh*k1/2);
    k2 = D2y_handle(x(end) + parameterh/2,...
        y(end) + parameterh*m1/2, Dy(end) + parameterh*k1/2);
    
    m3 = Dy_handle(x(end) + parameterh/2,...
        y(end) + parameterh*m2/2, Dy(end) + parameterh*k2/2);
    k3 = D2y_handle(x(end) + parameterh/2,...
        y(end) + parameterh*m2/2, Dy(end) + parameterh*k2/2);
    
    m4 = Dy_handle(x(end) + parameterh,...
        y(end) + parameterh*m3, Dy(end) + parameterh*k3);
    k4 = D2y_handle(x(end) + parameterh,...
        y(end) + parameterh*m3, Dy(end) + parameterh*k3);
    
    x(end + 1) = x(end) + parameterh;
    Dy(end + 1) = Dy(end) + (k1 + 2*k2 + 2*k3 + k4)*parameterh/6;
    y(end + 1) = y(end) + (m1 + 2*m2 + 2*m3 + m4)*parameterh/6;
end
[output] = [x;y;Dy];
end
