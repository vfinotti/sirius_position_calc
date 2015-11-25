function [y] = amean(x)

xsize = length(x);
y = zeros(1,xsize);

for i = 1:xsize
    if i == xsize
        y(i) = y(i-1); %(x(i)+0)/2;
    else
        y(i) = (x(i)+x(i+1))/2;
    end
end
    