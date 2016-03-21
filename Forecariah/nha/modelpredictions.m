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
c = clock;
fix(c)
red = [0,0.25,0.50];
Runforward = {'continue','yearly vector control','yearly vector control with active case-finding','yearly vector control with alternate active case-finding'}
for s = 2:4
    r = 0.50;
    load initconds
    X(:,4) = (1-r)*mkr;
    intervention = Runforward{s}
    Y = [];
    for i = 1:length(B)
        [t,y] = runHATintervention(X(i,:),Q(i,:),intervention,0);
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

    for i = 1:length(Y)
        for j = 1:length(Y{1})
            if j <= 360
                inc(j,i) = (Y{i}(j,8)+Y{i}(j,9))*360/Y{i}(j,6);
            else
                inc(j,i) = (Y{i}(j,15)-Y{i}(j-360,15))/Y{i}(j-360,6);
            end
        end
    end



    % Calculate total prevalences at each time point till 2030;
    Shai = (cellfun(@(y) (y(:,8)+y(:,9))'./sum(y(:,6:10)'),Y,'UniformOutput',false));
    TotPrev = cell2mat(Shai');
    VP = cell2mat((cellfun(@(y) y(:,4)'./sum(y(:,2:5)'), Y, ...
                           'UniformOutput',false))');
    LP = cell2mat((cellfun(@(y) y(:,13)'./sum(y(:,11:14)'), Y,'UniformOutput',false))');


    for i = 1:size(TotPrev,2)
        Tot = inc(i,:);%TotPrev(:,i);
        W=L./sum(L);
        E1(i) = length(Tot(1000*Tot<1)); % Eradication criteria
        E2(i) = length(Tot(10000*Tot<1)); % Elimination as public
                                          % health c   riteria
        E3(i) = length(Tot(100000*Tot<1));
        P1(i) = E1(i)/length(Tot);
        P2(i) = sum(L(10000*inc(i,:)<1)./sum(L))*E2(i)/length(Tot);
        P3(i) = E3(i)/length(Tot);
    end
     if r== 0
       hmm = s;
    else
        hmm = (s+r)*100
    end


    filename = sprintf('PredictionsA%d.mat',hmm);
    save(filename,'t','inc','TotPrev','VP','LP','E1','E2','E3','P1','P2','P3')
    clear all
    Runforward = {'continue','yearly vector control','yearly vector control with active case-finding','yearly vector control with alternate active case-finding'}
end


%% Use this part to plot the projections of each strategies

%  load Predictions1
% fig1 = figure('Position',[100,100,1000,500]); % Total prevalences and vector prevalences
% subplot(3,1,1)
% plot(t,TotPrev,'color',[0.5,0.5,0.5])
% title('HAT prevalences')
% ax = gca;
% ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
% ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
% subplot(3,1,2)
% plot(t,VP,'color',[0.5,0.5,0.5])
% title('Vector prevalences')
% ax = gca;
% ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];
% ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
% subplot(3,1,3)
% plot(t,LP,'color',[0.5,0.5,0.5])
% title('Animal prevalences')
% ax = gca;
% ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];
% ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}





%% Fig 3: Plot probabilities of elimination
fig3 = figure('Position',[100,100,1200,800]);
ColorSet = varycolor(4);
set(gca,'ColorOrder',ColorSet);
hold all;
for i = 1:4
    filename = sprintf('PredictionsA%d.mat',i)
    load(filename)
    plot(t,P2,'linewidth',1.5)
    %    hold on;
end
title('Probability of elimination as public health problem')
legend('No control','Annual vector control','Annual vector control and active case-finding','Annual vector control and biennial active case-finding ')
% ,'Vector control scale-up by 50% ','Vector control scale-up by 75% ','Case-finding scale-down to 35% ','Case-finding scale-down to 30% ','Vector control and case finding scale-down to 35% ', 'Vector control and case finding scale-down to 30% ')
legend('Location','southeast')
legend('boxoff')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
hold off
box('off')



ind = [2,3,4]
j = 1;
for i = ind
    filename = sprintf('PredictionsA%d.mat',i)
    load(filename);
    [a,b] = min(abs(t-7));
    ProbMat(j,1) = P2(b)
    j = j+1;
end
j = 1;
for i = ind
    filename = sprintf('PredictionsA%d25.mat',i)
    load(filename);
    [a,b] = min(abs(t-7));
    ProbMat(j,2) = P2(b)
    j = j+1;
end
j = 1;
for i = ind
    filename = sprintf('PredictionsA%d50.mat',i)
    load(filename);
    [a,b] = min(abs(t-7));
    ProbMat(j,3) = P2(b)
    j = j+1;
end


% fig3 = figure('Position',[100,100,1200,800]);
% ColorSet = varycolor(5);
% set(gca,'ColorOrder',ColorSet);
% hold all;
% for i = 1:5
%     filename = sprintf('Predictions%d.mat',i)
%     load(filename)
%     plot(t,P2,'linewidth',1.5)
%     %    hold on;
% end
% title('Prob. of Elim. as Public health problem')
% legend('No control','Annual vector control','Annual active case-finding','Annual vector control and active case-finding','Annual vector control and biennial active case-finding ')
% % ,'Vector control scale-up by 50% ','Vector control scale-up by 75% ','Case-finding scale-down to 35% ','Case-finding scale-down to 30% ','Vector control and case finding scale-down to 35% ', 'Vector control and case finding scale-down to 30% ')
% legend('Location','southeast')
% legend('boxoff')
% ax = gca;
% ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
% ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
% hold off
% box('off')


%% Fig 4: Plot barchart for prob. of elim. before 2020

% ind1 = [2,6,7]; ind2 = [3,8,9]; ind3 = [4,10,11];
% %ind = vertcat(ind1,ind2,ind3);
% ProbMat = zeros(3,3);
% j = 1;
% for i = ind1
%     filename = sprintf('Predictions%d.mat',i)
%     load(filename);
%     [a,b] = min(abs(t-8));
%     ProbMat(j,1) = P2(b);
%     j = j+1;
% end
% j = 1;
% for i = ind2
%     filename = sprintf('Predictions%d.mat',i)
%     load(filename);
%     [a,b] = min(abs(t-8));
%     ProbMat(j,2) = P2(b);
%     j = j+1;
% end
% j = 1;
% for i = ind3
%     filename = sprintf('Predictions%d.mat',i)
%     load(filename);
%     [a,b] = min(abs(t-8));
%     ProbMat(j,3) = P2(b);
%     j = j+1;
% end

% fig4 = figure('Position',[100,100,1200,800]);
% bar_handle = bar(ProbMat')
% set(bar_handle(1),'FaceColor',[252,141,89]/255);
% set(bar_handle(2),'FaceColor',[252,255,191]/255);
% set(bar_handle(3),'FaceColor',[145,191,219]/255);
% box('off')
% ylabel('Probability of HAT elimination by the end of 2020')
% ax = gca;
% %ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
% ax.XTickLabel = {'Annual vector Control','Annual case-finding','Annual vector control with case-finding ' }
