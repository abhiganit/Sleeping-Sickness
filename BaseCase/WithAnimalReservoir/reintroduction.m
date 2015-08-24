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
%s = 4;
for s = 4:5
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
        E(i) = length(Tot(10000*Tot>1)); % Elimination as public
                                          % health criteria
        P(i) = E(i)/length(Tot);
    end


    filename = sprintf('Reentry%d.mat',s);
    save(filename,'t','TotPrev','VP','E','P')
    clear all
    Runforward = {'continue','yearly vector control','only yearly active case-finding','yearly vector control with active case-finding','yearly vector control with alternate active case-finding','vector-control scale-up 1','vector-control scale-up 2','active case finding scale-down 1','active case finding scale-down 2','yearly vector control with active case-finding scale-up 1','yearly vector control with active case-finding scale-up 2'};
end


%% Ignore this part (combining the old data with new)
ColorSet = varycolor(2);
set(gca,'ColorOrder',ColorSet);
hold all;
for i = 1:2
    filename = sprintf('Reentry%d.mat',i+3);
    load(filename)
    plot(t,P);
    %    hold on;
end
hold off;


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
