function [output] = inv_DFT(x_values,max_frequency,stepsize)

% Given a set of data and max_frequency, inv_DFT evaluates the inverse 
% Discrete Fourier Transform of the given data set and has an output of 
% amplitudes X and frequencies k.
freq = 1:stepsize:max_frequency;
for n = 1:length(x_values)
    for k = 1:length(freq)
        f(k) = x_values(k)*exp(i*2*pi*(freq(k)-1)*(n-1)/(length(x_values)+1));
    end
    x(n) = sum(f)/length(x_values);
    t_values(n) = n;
end
output = [x;t_values];
end