
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


ProbMatA = zeros(3,3);
j = 1;
for i = ind1
    filename = sprintf('PredictionsA%d.mat',i)
    load(filename);
    [a,b] = min(abs(t-8));
    ProbMatA(j,1) = P2(b);
    j = j+1;
end
j = 1;
for i = ind2
    filename = sprintf('PredictionsA%d.mat',i)
    load(filename);
    [a,b] = min(abs(t-8));
    ProbMatA(j,2) = P2(b);
    j = j+1;
end
j = 1;
for i = ind3
    filename = sprintf('PredictionsA%d.mat',i)
    load(filename);
    [a,b] = min(abs(t-8));
    ProbMatA(j,3) = P2(b);
    j = j+1;
end


fig4 = figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
bar_handle = bar(ProbMat')
set(gca,'XTick',[])
set(bar_handle(1),'FaceColor',[252,141,89]/255);
set(bar_handle(2),'FaceColor',[252,255,191]/255);
set(bar_handle(3),'FaceColor',[145,191,219]/255);
box('off')
ylabel('Probability of HAT elimination by the end of 2020')
subplot(2,1,2)
bar_handle = bar(ProbMatA')
set(bar_handle(1),'FaceColor',[252,141,89]/255);
set(bar_handle(2),'FaceColor',[252,255,191]/255);
set(bar_handle(3),'FaceColor',[145,191,219]/255);
box('off')
ax = gca;
ax.XTickLabel = {'Annual vector Control','Annual case-finding','Annual vector control with case-finding' }
[hx,hy] = format_ticks(gca,{'Vector control',...
                    'Vector control  (50% scale-up)',...
                    'Vector control (75% scale-up)',...
                    'Case-finding',...
                    'Case-finding (35% coverage)',...
                    'Case-finding (30% coverage)',...
                    'Vector control with case-finding ',...
                    'Vector control with case-finding (35% coverage)',...
                    'Vector control with case-finding (30% coverage)'},...
                      [],[0.75,1,1.21,1.75,2,2.21,2.75,3,3.21],[], ...
                       45,0);

% ax = gca;
% ax.XTick = [0.75,1,1.21,1.75,2,2.21,2.75,3,3.21];
% ax.XTickLabel = {'Annual vector control',...
%                     'Vector control  (50% scale-up)',...
%                     'Vector control (75% scale-up)',...
%                     'Case-finding',...
%                     'Case-finding (35% coverage)',...
%                     'Case-finding (30% coverage)',...
%                     'Vector control with case-finding ',...
%                     'Vector control with case-finding (35% coverage)',...
%                     'Vector control with case-finding (30% coverage)'},...

% fix_xticklabels();