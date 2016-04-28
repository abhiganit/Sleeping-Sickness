function[Data] = BayesianMCMC(HealthZone,seed)
Data = xlsread('DRCdatasheet', HealthZone); % Year/Active coverage/New
rand('state',seed)
randn('state',seed)
load(HealthZone)

% Initialize
theta0 = ics(seed,:);
init_theta = theta0;
prior0 = calcprior(theta0);
out = runHATmodel(theta0,Data);
LogLik0 = out{4}
iterate  = 0;
scale = 0.05; %0.001*(2.38/sqrt(6));
Sigmaold = eye(6);
Sigmanew = Sigmaold;
Sigma = Sigmaold;
accept = 0;
for i = 1:2000
    %theta1 = mvnrnd(theta0,scale*Sigma);
    theta1 = theta0 - scale*init_theta.*(ones(1,6)-2*rand(1,6));
    prior1 = calcprior(theta1)
    out = runHATmodel(theta1,Data);
    LogLik1 = out{4};
    if (LogLik1 == -inf) || (prior1 == 0)
        i = i-1;
    else
        ratio = exp(theta1-theta0); %theta1/theta0;
        if rand < ratio
            accept = accept + 1
            theta0 = theta1;
            LogLik1 = LogLik0;
        end
        saveparams(i,:) = theta0;
        
    end
    save('posterior','saveparams')
end





