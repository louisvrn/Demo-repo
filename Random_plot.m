%% Clear workspace
clear 
close all
clc
%% Stupid plot 

for i = 1:1000
    x(i) = rand;
    y(i) = rand;
end
plot(x,y,'ro');