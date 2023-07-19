function [output] = DFT(x_values,max_frequency,stepsize)

% Given a set of data and max_frequency, DFT evaluates the Discrete Fourier
% Transform of the given data set and has an output of amplitudes X and
% frequencies k.

freq = 1:stepsize:max_frequency;
for k = 1:length(freq)
    for n = 1:length(x_values)
        f(n) = x_values(n)*exp(-i*2*pi*(freq(k)-1)*(n-1)/(length(x_values)+1));
    end
    X(k) = sum(f);
    k_values(k) = freq(k);
end
output = [X;k_values];
end