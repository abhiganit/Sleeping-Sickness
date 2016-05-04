function[] = Modelfitting13()
fix(clock)
[stupid, sheetname] = xlsfinfo('DRCdatasheet');
name = sheetname{13};
Data = xlsread('DRCdatasheet',name); % Year/Active coverage/New
                                         % cases(AC)/New cases(PC)/Population

 for i = 1:length(Data(:,3));
     X = betarnd(Data(i,3),ceil(Data(i,2).*Data(i,5)),10000,1)*ceil(Data(i,2).*Data(i,5));
      CIA(i,:) = quantile(X,[0,1]);
 end

N = 5000000  ;

parfor j = 1:N
    betaVH = 5*rand;
    nuH = 0.8*rand;
    k =  0.01*rand;
    betaH = 5*rand;
    zeta = 10*rand;
    a = Data(:,2);
    scal = rand;
    x0 = [betaVH,k,nuH,betaH,zeta,scal];
    param(j,:) = x0;
    out = runHATmodel(x0,Data,CIA);
    PD(j,:) = out{1};
    AD(j,:) = out{2};
    LogLik(j) = out{4};
end

%% Non-zero Likelihoods and respective parameters and outcomes
Ind = find(LogLik~=-Inf);
param1 = param(Ind,:);
Lik = exp(LogLik(Ind));
[val,ind] = max(Lik)

%[X,I] = sort(Lik,'descend')


fix(clock)
%sprintf('No of non-zero likelihoods = %d',length(Ind))
for i = 1:length(Ind)
    StP(i,:) = PD((Ind(i)),:);
    StA(i,:) = AD((Ind(i)),:);
    Pars(i,:) = param(Ind(i),:);
end

save(name,'val','ind','StP','StA','Pars')

end
