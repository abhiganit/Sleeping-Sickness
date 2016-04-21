%global Data
Data = xlsread('DRCdatasheet','Banga-Bola'); % Year/Active coverage/New
                                         % cases(AC)/New cases(PC)/Population
                                         
                                         
% x0 = [3.2910    0.0411    0.0088    4.6162    0.3105];
% out = runHATmodel1(x0,Data)
N = 50000;
parfor j = 1:N
    betaVH = 5*rand;
    nuH = 0.5*rand;
    k = 0.01*rand;
    betaH = 5*rand;
    zeta = 0.7*rand;
    a = Data(:,2);
    scal = rand;  %rand*(1/max(a(a~=0)));
    x0 = [betaVH,k,nuH,betaH,zeta,scal];
    param(j,:) = x0;
    Lik(j) = runHATmodel(x0,Data)
end

Ind = find(Lik~=-Inf)
NLIk = Lik(Ind)

for i = 1:length(Ind)
    out = runHATmodel1(param(Ind(i),:),Data);
    StP(i,:) = out{1};
    StA(i,:) = out{2};
end

PC = StP'
AC = StA'
t  = 1:length(Data(:,3));
fig1 = figure;
plot(t,AC(:,AC(1,:)<200),'color',[0.7,0.7,0.7])
hold on;
plot(t,Data(:,3),'ro')
hold off;

fig2 = figure;
plot(t,PC(:,AC(1,:)<200),'color',[0.7,0.7,0.7]);
hold on;
plot(t,Data(:,4),'go')
hold off;





%x = fminsearch(@runHATmodel,x0)


% Priors







% Loglikelihood






% MCMC steps