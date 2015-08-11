%% Ignore this part (generating data to plot before 2008 + after 2013 together)
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


%% Runforward loops over all the controls and saves data in Predictions%.mat

Runforward = {'continue','yearly vector control','only yearly active case-finding','yearly vector control with active case-finding','yearly vector control with alternate active case-finding','vector-control scale-up 1','vector-control scale-up 2','active case finding scale-down 1','active case finding scale-down 2','yearly vector control with active case-finding scale-up 1','yearly vector control with active case-finding scale-up 2'}

%s = 7;
%s = 3;
for s = 11:length(Runforward)
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
    Runforward = {'continue','yearly vector control','only yearly active case-finding','yearly vector control with active case-finding','yearly vector control with alternate active case-finding','vector-control scale-up 1','vector-control scale-up 2','active case finding scale-down 1','active case finding scale-down 2','yearly vector control with active case-finding scale-up 1','yearly vector control with active case-finding scale-up 2'}
end


%% Ignore this part (combining the old data with new)
 % t = t + 5*ones(size(t));
 % t = horzcat(ts,t');
 % TotPrev = horzcat(TotPrev0,TotPrev);
 % VP = horzcat(V0,VP);



%% Use this part to plot the projections of each strategies

% load Predictions1;
% fig1 = figure('Position',[100,100,1000,500]); % Total prevalences and vector prevalences
% subplot(2,1,1)
% plot(t,TotPrev,'color',[0.5,0.5,0.5])
% title('HAT prevalences')
% ax = gca;
% ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
% ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
% subplot(2,1,2)
% plot(t(),VP,'color',[0.5,0.5,0.5])
% title('Vector prevalences')
% ax = gca;
% ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];
% ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}





%% Fig 3: Plot probabilities of elimination
fig3 = figure('Position',[100,100,1200,800]);
ColorSet = varycolor(5);
set(gca,'ColorOrder',ColorSet);
hold all;
for i = 1:5
    filename = sprintf('Predictions%d.mat',i)
    load(filename)
    plot(t,P2,'linewidth',1.5)
    %    hold on;
end
title('Prob. of Elim. as Public health problem')
legend('No control','Annual vector control','Annual active case-finding','Annual vector control and active case-finding','Annual vector control and biennial active case-finding ')
% ,'Vector control scale-up by 50% ','Vector control scale-up by 75% ','Case-finding scale-down to 35% ','Case-finding scale-down to 30% ','Vector control and case finding scale-down to 35% ', 'Vector control and case finding scale-down to 30% ')
legend('Location','southeast')
legend('boxoff')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
hold off
box('off')


%% Fig 4: Plot barchart for prob. of elim. before 2020

ind1 = [2,6,7]; ind2 = [3,8,9]; ind3 = [4,10,11];
%ind = vertcat(ind1,ind2,ind3);
ProbMat = zeros(3,3);
j = 1;
for i = ind1
    filename = sprintf('Predictions%d.mat',i)
    load(filename);
    [a,b] = min(abs(t-8));
    ProbMat(j,1) = P2(b);
    j = j+1;
end
j = 1;
for i = ind2
    filename = sprintf('Predictions%d.mat',i)
    load(filename);
    [a,b] = min(abs(t-8));
    ProbMat(j,2) = P2(b);
    j = j+1;
end
j = 1;
for i = ind3
    filename = sprintf('Predictions%d.mat',i)
    load(filename);
    [a,b] = min(abs(t-8));
    ProbMat(j,3) = P2(b);
    j = j+1;
end

fig4 = figure('Position',[100,100,1200,800]);
bar_handle = bar(ProbMat')
set(bar_handle(1),'FaceColor',[252,141,89]/255);
set(bar_handle(2),'FaceColor',[252,255,191]/255);
set(bar_handle(3),'FaceColor',[145,191,219]/255);
box('off')
ylabel('Probability of HAT elimination by the end of 2020')
ax = gca;
%ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
ax.XTickLabel = {'Annual vector Control','Annual case-finding','Annual vector control with case-finding' }


t = linspace(1,3,2);
tint = linspace(3.1,5,2);

for j = 1:4
tl = tint;
t = horzcat(t,tl);
tint = linspace(tint(end)+0.1,tint(end)+2,2);
end
t

252,141,89
255,255,191
145,191,219
