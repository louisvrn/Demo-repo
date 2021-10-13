function [x,y] = RandomPlot(a)
validateattributes(a, {'numeric'}, {'size', [1,1]})
%Function that generates two random arrays of length a
for i = 1:a
    x(i) = rand;
    y(i) = rand;
    
end

