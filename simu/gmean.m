function [y] = gmean(x)

xsize = length(x);
y = zeros(1,xsize);

for i = 1:xsize
    if i == xsize
        y(i) = y(i-1);%sqrt(x(i)*1);
    else
        y(i) = sqrt(x(i)*x(i+1));
    end
end
    