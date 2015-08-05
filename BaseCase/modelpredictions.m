load initconds
intervention = 'continue'
% Note: B is set of posterior parameters and Q is set of initial conditions
T = [];
Y = [];

for i = 1:255
    [t,y] = runHATintervention(X(i,:),Q(i,:),intervention);
    T{i} = t;
    Y{i} = y;
end


for i =1:255
    y = Y{i};
    S1(i) = (y(end,8))/sum(y(end,6:end));
    S2(i) = (y(end,9))/sum(y(end,6:end));
    Tot(i) = (y(end,8)+y(end,9))/sum(y(end,6:end));
    V(i) = y(end,4)/sum(y(end,2:5));
end

% High intensity (more than 1 cases per 1000)
% Medium intensity (more than 1 cases per 10000)
% Low intensity (more than 1 cases per 1000000)

E1 = length(Tot(1000*Tot<1))
E2 = length(Tot(10000*Tot<1))
E3 = length(Tot(1000000*Tot<1))
P1 = E1/length(Tot)
P2 = E2/length(Tot)
P3 = E3/length(Tot)
