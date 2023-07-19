function dy = EXhoResonanceFunc(t,y,P)
% Damped driven HO
% d^2xdt^2 = -((P.wo)^2)*x - P.gamma*dx/dt + (P.A)*sin(P.w*t) 
% Note: y(1) = x, y(2)= dx/dt
dy = zeros(2,1);    % A column vector to be returned

dy(1) = y(2); 
dy(2) = -((P.wo)^2)*y(1)- P.gamma*y(2)+ (P.A)*sin(P.w*t);

