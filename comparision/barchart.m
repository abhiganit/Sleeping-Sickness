
%% Fig 4: Plot barchart for prob. of elim. before 2020
P1 = [1.0000    0.8228    0.0098
            1.0000    1.0000    0.6206
            1.0000    1.0000    0.1596]; % Updated base-case

P2 = [0.7693 0.6668 0.4910
      0.9832 0.8092 0.5340
      0.8943 0.7491 0.5340]; % Updated NHA

E1 = [0     0     0
      0     0     0
      0     0     0];

E2 = [0.0088         0         0
    0.3158    0.0439         0
    0.0965    0.0263         0];

# Close all figures
close all
fig4 = figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
bar_handle = bar(P1)
set(gca,'XTick',[])
set(bar_handle(1),'FaceColor',[252,141,89]/255);
set(bar_handle(2),'FaceColor',[252,255,191]/255);
set(bar_handle(3),'FaceColor',[145,191,219]/255);
leg = legend('0% reduction','25% reduction','50% reduction')
a = get(gca,'TickLabel')
set(gca,'TickLabel',a,'FontSize',12)
set(leg,'Location','northeast');%'boxoff',
set(leg,'FontSize',14)
legend('boxoff')

box('off')
ylabel('Probability of HAT elimination by the end of 2020','FontSize',16)
subplot(2,1,2)
bar_handle = bar(P2)
set(bar_handle(1),'FaceColor',[252,141,89]/255);
set(bar_handle(2),'FaceColor',[252,255,191]/255);
set(bar_handle(3),'FaceColor',[145,191,219]/255);
box('off')
ax = gca;
ax.XTickLabel = {'Annual vector control',...
                 'Annual vector control with annual mass screening  ',...
                 'Annual vector control with biennial mass screening' }
a = get(gca,'TickLabel')
set(gca,'TickLabel',a,'FontSize',12)
fix_xticklabels(gca,0.1,{'FontSize',14});

close all
fig5 = figure('units','normalized','outerposition',[0 0 1 0.8])
A = [2023,2025,2030;2020,2022,2025;2021,2023,2027];
% A(:,3) = A(:,3)-A(:,2);
% A(:,2) = A(:,2)-A(:,1);
bar_handle = bar(A,'grouped')
set(gca,'XTick',[])
set(bar_handle(1),'FaceColor',[252,141,89]/255);
set(bar_handle(2),'FaceColor',[252,255,191]/255);
set(bar_handle(3),'FaceColor',[145,191,219]/255);
box('off')
ylim([2013,2030])
leg = legend(' 0% reduction','25% reduction','50% reduction')
% ,'Vector control scale-up by 50% ','Vector control scale-up by 75% ','Case-finding scale-down to 35% ','Case-finding scale-down to 30% ','Vector control and case finding scale-down to 35% ', 'Vector control and case finding scale-down to 30% ')
set(leg,'Location','northeast');%'boxoff',
set(leg,'FontSize',14)
legend('boxoff')
ylabel('Years','FontSize',16)
ax = gca;
ax.XTick = [1,2,3]
ax.XTickLabel = {'Annual vector control',...
                 'Annual vector control with annual case-finding ',...
                 'Annual vector control with biennial case-finding' }
a = get(gca,'TickLabel')
set(gca,'TickLabel',a,'FontSize',12)
fix_xticklabels(gca,0.1,{'FontSize',14});
% [hx,hy] = format_ticks(gca,{'Vector control',...
%                     'Vector control  (50% scale-up)',...
%                     'Vector control (75% scale-up)',...
%                     'Case-finding',...
%                     'Case-finding (35% coverage)',...
%                     'Case-finding (30% coverage)',...
%                     'Vector control with case-finding ',...
%                     'Vector control with case-finding (35% coverage)',...
%                     'Vector control with case-finding (30% coverage)'},...
%                       [],[0.75,1,1.21,1.75,2,2.21,2.75,3,3.21],[], ...
%                        45,0);



%%% Not Using these right now;
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


% ProbMatA = zeros(3,3);
% j = 1;
% for i = ind1
%     filename = sprintf('PredictionsA%d.mat',i)
%     load(filename);
%     [a,b] = min(abs(t-8));
%     ProbMatA(j,1) = P2(b);
%     j = j+1;
% end
% j = 1;
% for i = ind2
%     filename = sprintf('PredictionsA%d.mat',i)
%     load(filename);
%     [a,b] = min(abs(t-8));
%     ProbMatA(j,2) = P2(b);
%     j = j+1;
% end
% j = 1;
% for i = ind3
%     filename = sprintf('PredictionsA%d.mat',i)
%     load(filename);
%     [a,b] = min(abs(t-8));
%     ProbMatA(j,3) = P2(b);
%     j = j+1;
% end
