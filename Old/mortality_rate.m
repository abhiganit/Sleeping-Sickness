function[out] = mortality_rate(rho, t,l,d)

l = l/12; % time in year
d = d/12; % time in year
t = mod(t,1); % bringing it between 0 and 1.

const_rate = rho;

if t <= l
    out = const_rate;
elseif t > l & t < l+d
    out = (const_rate./d)*(l+d-t);
else
    out = 0;
end
