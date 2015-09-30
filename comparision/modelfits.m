

fig1 = figure('units','normalized','outerposition',[0 0 0.7 1]) %figure();
load plotbasecase
t = [1,3,5,6]; % [2008,2010,2012,2013]
X = [t flip(t)];
X1 = [t flip(t)];
subplot(2,2,1)
%plot(t,B(:,1:4),'Color',[0.75,0.75,0.75]);
fill(X,Y,[0.75,0.75,0.75],'EdgeColor','none')
hold on;
plot(t,Best(1:4),'r','linewidth',1);
errorbar(t,DS1, bd1(:,1),bd1(:,2),'ko','linewidth',1.5)
ylabel('Prevalence','FontSize',16)
xlabel('Years','FontSize',16)
%plot(tdv,DV1,'ko','linewidth',1.5);
hold off;
title('Stage I','FontSize',16)
ax = gca;
ax.XTick = [1,2,3,4,5,6]
ax.XTickLabel = {'2008','2009', '2010','2011','2012',' 2013'}
a = get(gca,'TickLabel')
set(gca,'TickLabel',a,'FontSize',12)

box('off')
subplot(2,2,2)
%plot(t,B(:,5:8),'Color',[0.75,0.75,0.75]);
fill(X1,Y1,[0.75,0.75,0.75],'EdgeColor','none')
hold on;
plot(t,Best(5:8),'r','linewidth',1);
errorbar(t,DS2, bd2(:,1),bd2(:,2),'ko','linewidth',1.5)
%plot(tdv,DV2,'ko','linewidth',1.5);
xlabel('Years','FontSize',16)
hold off;
title('Stage II','FontSize',16)
ax = gca;
ax.XTick = [1,2,3,4,5,6]
ax.XTickLabel = {'2008','2009', '2010','2011','2012',' 2013'}
a = get(gca,'TickLabel')
set(gca,'TickLabel',a,'FontSize',12)


box('off')
clear all
load plotNHA
t = [1,3,5,6]; % [2008,2010,2012,2013]
X = [t flip(t)];
X1 = [t flip(t)];
subplot(2,2,3)
%plot(t,B(:,1:4),'Color',[0.75,0.75,0.75]);
fill(X,Y,[0.75,0.75,0.75],'EdgeColor','none')
hold on;
plot(t,Best(1:4),'r','linewidth',1);
errorbar(t,DS1, bd1(:,1),bd1(:,2),'ko','linewidth',1.5)
ylabel('Prevalence','FontSize',16)
xlabel('Years','FontSize',16)
%plot(tdv,DV1,'ko','linewidth',1.5);
hold off;
title('Stage I','FontSize',16)
ax = gca;
ax.XTick = [1,2,3,4,5,6]
ax.XTickLabel = {'2008','2009', '2010','2011','2012',' 2013'}
a = get(gca,'TickLabel')
set(gca,'TickLabel',a,'FontSize',12)
box('off')
subplot(2,2,4)
%plot(t,B(:,5:8),'Color',[0.75,0.75,0.75]);
fill(X1,Y1,[0.75,0.75,0.75],'EdgeColor','none')
hold on;
plot(t,Best(5:8),'r','linewidth',1);
errorbar(t,DS2, bd2(:,1),bd2(:,2),'ko','linewidth',1.5)
%plot(tdv,DV2,'ko','linewidth',1.5);
xlabel('Years','FontSize',16)
hold off;
title('Stage II','FontSize',16)
ax = gca;
ax.XTick = [1,2,3,4,5,6]
ax.XTickLabel = {'2008','2009', '2010','2011','2012',' 2013'}
a = get(gca,'XTickLabel')
set(gca,'XTickLabel',a,'FontSize',12)
%fix_xticklabels(gca,0.1,{'FontSize',14});
box('off')
