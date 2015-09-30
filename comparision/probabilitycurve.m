
close all
fig3 = figure('units','normalized','outerposition',[0 0 1 1]) %figure('Position',[100,100,1200,800]);
P = [27,158,119
217,95,2
117,112,179
231,41,138]/255;
subplot(2,1,1)
% ColorSet = varycolor(4);
% set(gca,'ColorOrder',ColorSet);
hold all;
for i = 1:4
    filename = sprintf('Predictions%d.mat',i)
    load(filename)
    plot(t,P2,'linewidth',1.5,'Color',P(i,:))
      hold on;
end
%xlabel('Years','FontSize',16)
ylabel('Probability of elimination as public health problem','FontSize',16)
leg = legend('No control','Annual vector control','Annual vector control with annual mass screening','Annual vector control with biennial mass screening ')
% ,'Vector control scale-up by 50% ','Vector control scale-up by 75% ','Case-finding scale-down to 35% ','Case-finding scale-down to 30% ','Vector control and case finding scale-down to 35% ', 'Vector control and case finding scale-down to 30% ')
set(leg,'Location','southeast');%'boxoff',
set(leg,'FontSize',14)
legend('boxoff')
% legend('Location','southeast')
% legend('boxoff')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019', ...
                 '2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
a = get(gca,'TickLabel')
set(gca,'TickLabel',a,'FontSize',12)
hold off
box('off')
subplot(2,1,2)
% ColorSet = varycolor(4);
% set(gca,'ColorOrder',ColorSet);
hold all;
for i = 1:4
    filename = sprintf('PredictionsA%d.mat',i)
    load(filename)
    plot(t,P2,'linewidth',1.5,'Color',P(i,:))
      hold on;
end
xlabel('Years','FontSize',16)
%ylabel('Probability of elimination as public health problem','FontSize',16)
% legend('No control','Annual vector control','Annual vector control with annual mass screening','Annual vector control with biennial mass screening ')
% % ,'Vector control scale-up by 50% ','Vector control scale-up by 75% ','Case-finding scale-down to 35% ','Case-finding scale-down to 30% ','Vector control and case finding scale-down to 35% ', 'Vector control and case finding scale-down to 30% ')
% legend('Location','northeast')
% legend('boxoff')
ax = gca;
ax.XTick = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
ax.XTickLabel = {'2013','2014','2015','2016','2017','2018','2019', ...
                 '2020','2021','2022','2023','2024','2025','2026','2027','2028','2029','2030'}
a = get(gca,'TickLabel')
set(gca,'TickLabel',a,'FontSize',12)
hold off
box('off')
