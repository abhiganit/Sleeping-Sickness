function[Data] = BayesianMCMC(HealthZone,seed)
Data = xlsread('DRCdatasheet', HealthZone); % Year/Active coverage/New
rand('state',seed)
randn('state',seed)
load(HealthZone)

for i = 1:length(Data(:,3));
     X = betarnd(Data(i,3),ceil(Data(i,2).*Data(i,5)),10000,1)*ceil(Data(i,2).*Data(i,5));
      CIA(i,:) = quantile(X,[0,1]);
end

% Initialize
theta0 = Pars(ind,:); %ics(seed,:);
init_theta = theta0;
prior0 = calcprior(theta0);
out = runHATmodel(theta0,Data,CIA);
LogLik0 = out{4}
scale = 0.001*(2.38/sqrt(6));
Sigmaold = eye(6);
Sigmanew = Sigmaold;
Sigma = Sigmaold;
accept = 0;
sizeofslot = 500;
total = 500000
iterate = 1;
while (iterate<total)
    accept = 0;
    slot = 1;
    Sigma = scale*(0.25*Sigmaold+0.75*Sigmanew);
    while slot < sizeofslot
        theta1 = mvnrnd(theta0,scale*Sigma);
        %theta1 = theta0 - scale*init_theta.*(ones(1,6)-2*rand(1,6));
        prior1 = calcprior(theta1);
        out = runHATmodel(theta1,Data,CIA);
        LogLik1 = out{4};
        if (LogLik1 ~= -inf)
            ratio = exp(theta1-theta0); %theta1/theta0;
            if rand < ratio
                accept = accept + 1;
                theta0 = theta1;
                LogLik1 = LogLik0;
            end
            thetaslot(slot,:) = theta0;
            saveparams(iterate,:) = theta0;
        end
        iterate = iterate+1;
        slot = slot+1;
    end
    accept_prob = accept/sizeofslot
    Sigmaold = Sigmanew;
    Sigmanew = cov(thetaslot);
    scale = scale*exp(accept_prob-0.23)/(iterate+1);
    scale';

    save('posterior','saveparams')
end
