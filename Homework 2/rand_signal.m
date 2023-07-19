function [output] = rand_signal(t_range,t_stepsize,max_amplitude,noise)

% Function rand_signal outputs a sum of different sinusoids and added
% noise.

total_sinusoid = randi([4 15],1);
t = t_range(1):t_stepsize:t_range(2);
for count = 1:total_sinusoid
    sinusoid(count,:) = max_amplitude*rand(1)*cos((rand(1)+1)*t ...
        + (rand(1)*2*pi));
end
if noise == 1
    for count = 1:length(t)
        n(count) = randn(1)*1.5;
    end
    sinusoid(end+1,:) = n;
end
x_data = sum(sinusoid);
t_data = t;
[output] = [t_data;x_data];