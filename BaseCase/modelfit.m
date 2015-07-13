function[posterior,params] = modelfit()

%% Bayesian Melding
    tic;
    N = 100000;
    % N samples from priors
    params = zeros(N,3);

    parfor i = 1:N
        betaVH = 10*rand;
        betaH = 10*rand;
        zeta2 = 2.28*rand;
        params(i,:)  = [betaVH,betaH,zeta2];

        [a,b,c] = runHATmodel(params(i,:));

        Likelihood(i) = betapdf(a,5,4307)*betapdf(b,7,4307)*betapdf(c,1,874);
    end
    toc

    parfor i = 1:N
    weights(i) = Likelihood(i)/sum(Likelihood);

    end

    weight = weight(weight~=0);
    params = params(weight~=0,:);

    % {nter1 = Desired number of sample for poserior
       k1 = 1 ;

    %}

    total = 500;
    j = 1;
    while j < total
     k =randi(length(weight),1);

     if(rand <weight(k))
     Par =params(k,:);
     posterior(j,:)=Par;
     j = j+1
    end
end
