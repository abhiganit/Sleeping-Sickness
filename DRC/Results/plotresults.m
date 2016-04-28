
[~,healthzones] = xlsfinfo('DRCdatasheet');

filename = healthzones{25}
Data = xlsread('DRCdatasheet',filename);
load(filename)



PC = StP';
AC = StA';
t  = 1:length(Data(:,3));
fig1 = figure;
plot(t,AC(:,AC(1,:)<400),'color',[0.7,0.7,0.7])
hold on;
plot(t,Data(:,3),'ro',t,AC(:,ind),'k','linewidth',2)
title('Active Cases')
hold off;

fig2 = figure;
plot(t,PC(:,AC(1,:)<400),'color',[0.7,0.7,0.7]);
hold on;
plot(t,Data(:,4),'go',t,PC(:,ind),'k','linewidth',2)
title('Passive Cases')
hold off;
