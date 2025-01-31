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
Runforward = {'continue','yearly vector control','yearly vector control with active case-finding','yearly vector control with alternate active case-finding'};
for s = 2:4
    r = 0.50; % Percentage reduction in maximal kill rate.
    load initconds
    X(:,4) = (1-r)*mkr;
    intervention = Runforward{s}
    Y = [];
    for i = 1:length(B)
        [t,y] = runHATintervention(X(i,:),Q(i,:),intervention,r);
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
                inc(j,i) = (Y{i}(j,11)-Y{i}(j-360,11))/Y{i}(j-360,6);
            end
        end
    end


    % Calculate total prevalences at each time point till 2030;
    Shai = (cellfun(@(y) (y(:,8)+y(:,9)+y(:,10))'./sum(y(:,6:end)'),Y,'UniformOutput',false));
    TotPrev = cell2mat(Shai');




    VP = cell2mat((cellfun(@(y) y(:,4)'./sum(y(:,2:5)'), Y,'UniformOutput',false))');


    for i = 1:size(TotPrev,2)
        Tot = inc(i,:);
        E1(i) = length(Tot(1000000*Tot<1)); % Eradication criteria
        E2(i) = length(Tot(10000*Tot<1)); % Elimination as public
                                          % health c   riteria
        E3(i) = length(Tot(100000*Tot<1));
        P1(i) = sum(L(1000000*inc(i,:)<1)./sum(L))*E1(i)/length(Tot);
        P2(i) = sum(L(10000*inc(i,:)<1)./sum(L))*E2(i)/length(Tot);
        P3(i) = E3(i)/length(Tot);
    end
    if r== 0
       hmm = s;
    else
        hmm = (s+r)*100
    end
    filename = sprintf('Predictions%d.mat',hmm);
    save(filename,'t','inc','TotPrev','VP','E1','E2','E3','P1','P2','P3')
    clear all
    Runforward = {'continue','yearly vector control','yearly vector control with active case-finding','yearly vector control with alternate active case-finding'}
end




%% Use this part to plot the projections of each strategies

% load Predictions2;
% fig1 = figure('Position',[100,100,1000,500]); % Total prevalences and vector prevalences
% subplot(2,1,1)
% plot(t,inc,'color',[0.5,0.5,0.5])
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



% P = [166,97,26
% 223,194,125
% 128,205,193
% 1,133,113]/255

%% Fig 3: Plot probabilities of elimination
fig3 = figure('units','normalized','outerposition',[0 0 1 1]) %figure('Position',[100,100,1200,800]);
P = [230,97,1
253,184,99
178,171,210
94,60,153]/255;
ColorSet = varycolor(4);
set(gca,'ColorOrder',ColorSet);
hold all;
for i = 1:4
    filename = sprintf('Predictions%d25.mat',i)
    load(filename)
    plot(t,P3,'linewidth',1.5,'Color',P(i,:))
      hold on;
end
xlabel('Years','FontSize',16)
ylabel('Probability of elimination as public health problem','FontSize',16)
legend('No control','Annual vector control','Annual vector control with annual mass screening','Annual vector control with biennial mass screening ')
% ,'Vector control scale-up by 50% ','Vector control scale-up by 75% ','Case-finding scale-down to 35% ','Case-finding scale-down to 30% ','Vector control and case finding scale-down to 35% ', 'Vector control and case finding scale-down to 30% ')
legend('Location','northeast')
legend('boxoff')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
hold off
box('off')


%% Fig 4: Plot barchart for prob. of elim. before 2020


% ind1 = [2,6,7]; ind2 = [3,8,9]; ind3 = [4,10,11];
% %ind = vertcat(ind1,ind2,ind3);

ind = [2,3,4]
j = 1;
for i = ind
    filename = sprintf('Predictions%d.mat',i)
    load(filename);
    [a,b] = min(abs(t-7));
    ProbMat(j,1) = P2(b)
    j = j+1;
end
j = 1;
for i = ind
    filename = sprintf('Predictions%d25.mat',i)
    load(filename);
    [a,b] = min(abs(t-7));
    ProbMat(j,2) = P2(b)
    j = j+1;
end
j = 1;
for i = ind
    filename = sprintf('Predictions%d50.mat',i)
    load(filename);
    [a,b] = min(abs(t-7));
    ProbMat(j,3) = P2(b)
    j = j+1;
end
