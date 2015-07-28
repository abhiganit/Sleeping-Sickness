function[out] = coverage(t)

if t > 1 & t <= 2
    out = 0.68;
elseif t > 3 & t <=4  % 2012
    out = 0.46;
elseif t > 4 & t <= 5 % 2013
    out = 0.54;
else
    out = 0;
end
