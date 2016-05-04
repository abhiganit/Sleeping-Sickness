close all;
clear all;
[~,healthzones] = xlsfinfo('DRCdatasheet');

filename = healthzones{4}
% 3
Data = xlsread('DRCdatasheet',filename);
load(filename)

%%  Modelfit
% PC = StP';
% AC = StA';
% t  = 1:length(Data(:,3));
% fig1 = figure;
% plot(t,AC(:,AC(1,:)<400),'color',[0.7,0.7,0.7])
% hold on;
% plot(t,Data(:,3),'ro',t,AC(:,ind),'k','linewidth',2)
% title('Active Cases')
% hold off;

% fig2 = figure;
% plot(t,PC(:,AC(1,:)<400),'color',[0.7,0.7,0.7]);
% hold on;
% plot(t,Data(:,4),'go',t,PC(:,ind),'k','linewidth',2)
% title('Passive Cases')
% hold off;



for Strategy = 1:2;  % Maximum, Mean,
    parfor j = 1:length(Pars)
        out = runHATmodel(Pars(j,:),Data,Strategy);
        PD(j,:) = out{1};
        AD(j,:) = out{2};
        LogLik(j) = out{3};
    end
    Pas{Strategy} = PD;
    Act{Strategy} = AD;
    lglk{Strategy} = LogLik;
end


%% Bayesian melding step (Weighing and re-sampling)
Strategy = 1;
Lik = exp(lglk{Strategy})
parfor i = 1:length(Lik)
    weights(i) = Lik(i)/sum(Lik);
end

Total = 10000;
j = 1;
while j < Total+1
    k = randi(length(weights),1);
    if (rand<weights(k))
        posterior(j,:) = Pars(k,:);
        PD(j,:) = Pas{Strategy}(k,:);
        AD(j,:) = Act{Strategy}(k,:);

        j = j+1;
    end
end



PC = PD';
AC = AD';

Casesper10k = (AC+PC)*(10000/Data(1,5));

ProbofElim = sum((Casesper10k<1),2)/length(AC);

Last10{Strategy} = ProbofElim(end-10:end)


%% Results
t = 1:length(Data(:,2));
t1  = 1:length(Data(:,2))+18;
fig3 = figure;
plot(t1,AC(:,AC(1,:)<400),'color',[0.7,0.7,0.7])
hold on;
plot(t,Data(:,3),'ro',t1,AC(:,ind),'k','linewidth',2)
title('Active Cases')
hold off;

fig4 = figure;
plot(t1,PC(:,AC(1,:)<400),'color',[0.7,0.7,0.7]);
hold on;
plot(t,Data(:,4),'go',t1,PC(:,ind),'k','linewidth',2)
title('Passive Cases')
hold off;


% for Strategy = 1:2
%     PC = Pas{Strategy}';
%     AC = Act{Strategy}';

%     Casesper10k = (AC+PC)*(10000/Data(1,5));

%     ProbofElim = sum((Casesper10k<1),2)/length(AC);

%     Last10{Strategy} = ProbofElim(end-10:end);
% end

% Last10{1}


% Last10{2}
