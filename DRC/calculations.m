load posterior
params = saveparams;
name = 'Banga-Bola'
Data = xlsread('DRCdatasheet',name); % Year/Active coverage/New

for j = 1:length(params)
    out = runHATmodel(params(j,:),Data);
    PD(j,:) = out{1};
    AD(j,:) = out{2};
    LogLik(j) = out{4};
end

AC = AD'; PC = PD';
t  = 1:length(Data(:,3));
fig1 = figure;
plot(t,AC,'color',[0.7,0.7,0.7])
hold on;
plot(t,Data(:,3),'ro')
hold off;

fig2 = figure;
plot(t,PC,'color',[0.7,0.7,0.7]);
hold on;
plot(t,Data(:,4),'go')
hold off;