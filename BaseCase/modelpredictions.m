load initconds
intervention = 'continue'
% Note: X is set of posterior parameters
%       Q is set of initial conditions
%       B is set of 2012 prevalences

Y = [];
for i = 1:255
    [t,y] = runHATintervention(X(i,:),Q(i,:),intervention);
     Y{i} = y;
end

Fail = [];
for i = 1:255
    if size(Y{i},1)<length(Y{1})
       Fail = [Fail,i];
    end
end
Y(Fail) = [];
X(Fail,:) = [];
Q(Fail,:) = [];
B(Fail,:) = [];


% Calculate total prevalences at each time point till 2030;
Shai = (cellfun(@(y) (y(:,8)+y(:,9))'./sum(y(:,6:end)'),Y,'UniformOutput',false));
TotPrev = cell2mat(Shai');
VP = cell2mat((cellfun(@(y) y(:,4)'./sum(y(:,2:5)'), Y,'UniformOutput',false))');


for i = 1:size(TotPrev,2)
    Tot = TotPrev(:,i);
    E1(i) = length(Tot(1000*Tot<1)); % Eradication criteria
    E2(i) = length(Tot(10000*Tot<1)); % Elimination as public health criteria
    E3(i) = length(Tot(1000000*Tot<1));
    P1(i) = E1(i)/length(Tot);
    P2(i) = E2(i)/length(Tot);
    P3(i) = E3(i)/length(Tot);
end


% Plots

fig1 = figure('Position',[100,100,1000,500]); % Total prevalences and vector prevalences
                %subplot(1,2,1)
plot(t,TotPrev,'color',[0.7,0.7,0.7])
title('HAT prevalences')
ax = gca;
ax.XTick = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}


%plot(t,VP)




fig2 = figure('Position',[100,100,1000,500]);; % Continuous Probability of elimination as Public
               % Health problem.
plot(t,P2)
title('Prob. of Elim. as PHP')
ax = gca;
ax.XTick = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
