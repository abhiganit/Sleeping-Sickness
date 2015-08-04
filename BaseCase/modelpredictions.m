load initconds
intervention = 'continue'
% Note: B is set of posterior parameters and Q is set of initial conditions
T = [];
Y = [];

for i = 1:255
    [t,y] = runHATintervention(B(i,:),Q(i,:),intervention);
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


E1 = 0; E2 = 0; E3 = 0;
for i = 1:255
    T1 = 1000*Tot(i);
    T2 = 10000*Tot(i);
    T3 = 1000000*Tot(i);
    if T1 < 1
        E1 = E1+1; E2 = E2+1; E3 = E3+1;
    elseif T2 < 1
        E2 = E2+1; E3 = E3+1;
    elseif T3 < 1
        E3 = E3+1;
    end
end
