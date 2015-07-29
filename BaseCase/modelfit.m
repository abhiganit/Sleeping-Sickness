function[] = modelfit()
%% Sampling.

    tic;

    % Data & Sample: Years(2008/10/12/13). Columns:  Stage 1 | Stage 2
    % | Vector Prevalence

    Data = [3,2,9.53;
            4,8, 0  ;
            7,13,0  ;
            3,4, 0 ];
    SampSize = [1488,1488,1634;
                4514,4514,0;
                7708,7708,0;
                7788,7788,0];

    X = betarnd(3,1488,10000,1);
    Y = betarnd(2,1488,10000,1);
    ci1 = quantile(X,[0.025,0.975])
    ci2 = quantile(Y,[0.025,0.975])



    N = 100000;
    % N samples from priors
    params = zeros(N,4);
    tic
    parfor j = 1:N
        betaVH = 0.1+0.5*rand;
        betaH = rand;
        zeta = 1.37*rand;
        rho = 365*0.15*rand;

        params(j,:)  = [betaVH,betaH,zeta,rho];

        out = runHATmodel(params(j,:));

        if (out(1)<=ci1(1)) || (out(1)>=ci1(2)) || (out(5)<=ci2(1)) || (out(5)>=ci2(2)) || out(9) > 0.010
	  Likelihood(j) = 0;
        else
            Lik1 = 1; Lik2 = 1;
            for i = 1:4
                Lik1 = betapdf(out(i),Data(i,1),SampSize(i,1));
                Lik2 = betapdf(out(i+4),Data(i,2),SampSize(i,2));
            end
            Likelihood(j) = Lik1*Lik2*betapdf(out(9),Data(1,3),SampSize(1,3));
        end
    end
    toc
    save('output','params','Likelihood')

end
