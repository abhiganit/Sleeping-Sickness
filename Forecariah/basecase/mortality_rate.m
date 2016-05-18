function[out] = mortality_rate(rho, t,l,d)

l = l/12; % time in year
d = d/12; % time in year
t = mod(t,1); % bringing it between 0 and 1.

const_rate = rho;

% out = 1-(1-const_rate)*t;

if t <= l
    out = const_rate;
elseif t > l & t < l+d
    out = (const_rate./d)*(l+d-t);
else
    out = 0;
end
%mortality_rate(rho, t, l, d)



% timepoints = linspace(0,1,365);

% j = 1;
% for t  = timepoints
% if t <= l
%     out(j) = const_rate;
% elseif t > l & t < l+d
%     out(j) = (const_rate./d)*(l+d-t);
% else
%     out(j) = 0;
% end
% j = j+1;
% end


% j = 1
% alpha = 0.3
% for i = t
%     out(j) = 1-(1-alpha)*i;
%     j =  j+1
% end
