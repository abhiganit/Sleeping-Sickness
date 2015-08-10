% Note: X is set of posterior parameters
%       Q is set of initial conditions
%       B is set of 2012 prevalences
% Ys = [];
% for i = 1:length(B)
%     out = runHATmodel(X(i,:));
%      ts = out{3};
%      Ys{i} = out{4};
%      Vs{i} = out{5};
% end
% TotPrev0 = cell2mat(Ys');
% V0 = cell2mat(Vs');

%save('Precontroldata','TotPrev0','ts','V0')
%load Precontroldata;

Runforward = {'continue','yearly vector control','only yearly active case-finding','yearly vector control with active case-finding','yearly vector control with alternate active case-finding'}

%s = 5;
for s = 1:length(Runforward)
    load initconds
    intervention = Runforward{s}
    Y = [];
    for i = 1:length(B)
        [t,y] = runHATintervention(X(i,:),Q(i,:),intervention);
        Y{i} = y;
    end

    Fail = [];
    for i = 1:length(B)
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
        E2(i) = length(Tot(10000*Tot<1)); % Elimination as public
                                          % health c   riteria
        E3(i) = length(Tot(1000000*Tot<1));
        P1(i) = E1(i)/length(Tot);
        P2(i) = E2(i)/length(Tot);
        P3(i) = E3(i)/length(Tot);
    end


    filename = sprintf('Predictions%d.mat',s);
    save(filename,'t','TotPrev','VP','E1','E2','E3','P1','P2','P3')
    clear all
    Runforward = {'continue','yearly vector control','only yearly active case-finding','yearly vector control with active case-finding','yearly vector control with alternate active case-finding'}
end

% Plots

% Combining the old data with new
 % t = t + 5*ones(size(t));
 % t = horzcat(ts,t');
 % TotPrev = horzcat(TotPrev0,TotPrev);
 % VP = horzcat(V0,VP);

%load Predictions4

fig1 = figure('Position',[100,100,1000,500]); % Total prevalences and vector prevalences
subplot(2,1,1)
plot(t,TotPrev,'color',[0.5,0.5,0.5])
title('HAT prevalences')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
subplot(2,1,2)
plot(t(),VP,'color',[0.5,0.5,0.5])
title('Vector prevalences')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}


for i = 1:5
    filename = sprintf('Predictions%d.mat',i)
    load(filename)
    plot(t,P2,'linewidth',1.5)
    hold on;
end
title('Prob. of Elim. as PHP')
legend('No control','Annual vector control','Annual active case-finding','Annual vector control and active case-finding','Annual vector control and biennial active case-finding')
legend('Location','southeast')
legend('boxoff')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
hold off


fig2 = figure('Position',[100,100,1000,500]);; % Continuous Probability of elimination as Public
               % Health problem.
plot(t,P2)
title('Prob. of Elim. as PHP')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
ax.XTickLabel = {'2012','2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
