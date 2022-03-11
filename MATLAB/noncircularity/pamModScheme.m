%  These can lead to non-circular signals
% PAM scheme, having an initial phase of pi/4
% M = 8
data = randi([0 M-1],100,1);
modData = pammod(data,M,pi/8);
scatterplot(modData); grid on;