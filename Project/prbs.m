function [r] = prbs(N,c,p)
%PRBS Summary of this function goes here
r = zeros(N,1);
x = rand;
if x <= p
    y = c;
else
    y = -c;
end
r(1) = y;
for i = 2:N
    x = rand;
    if x <= p
    r(i) = r(i-1);
    else
    r(i) = -r(i-1);
    end
end
end

