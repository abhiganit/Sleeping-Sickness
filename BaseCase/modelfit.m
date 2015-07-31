function[] = modelfit()
%% Sampling.

    tic;

    % Data & Sample: Years(2008/10/12/13). Columns:  Stage 1 | Stage 2
    % Vector Prevalence

    % Data = [7,21,9.53  ;
    %         3,7, 0 ];
    % SampSize = [21,7708,1634;
    %             7,7788,0];

    Data = [7,13,9.53;3,4,0];
    SampSize = [7708,7708,1634;7788,7788,0];


    X = betarnd(7,7708,10000,1);
    Y = betarnd(13,7708,10000,1);
    ci1 = quantile(X,[0.025,0.975]) % S1, 2008
    ci2 = quantile(Y,[0.025,0.975]) % S2, 2008


    N = 100000;
    % N samples from priors
    params = zeros(N,4);
    tic
    parfor j = 1:N
        betaVH = 0.1+0.5*rand;
        betaH = rand;
        zeta = 1.37*rand;
        %        cov1 = rand;
        rho = 365*0.25*rand;

        params(j,:)  = [betaVH,betaH,zeta,rho];

        out = runHATmodel(params(j,:));

        if (out(1)<=ci1(1)) || (out(1)>=ci1(2)) || (out(3)<=ci2(1)) ...
                    || (out(3)>=ci2(2)) || out(5) > 0.010
            Likelihood(j) = 0;
        else
            Lik1 = 1; Lik2 = 1;
            for i = 1:2
                Lik1 = betapdf(out(i),Data(i,1),SampSize(i,1));
                Lik2 = betapdf(out(i+2),Data(i,2),SampSize(i,2));
            end
            Likelihood(j) = Lik1*Lik2*betapdf(out(5),Data(1,3),SampSize(1,3));
        end
    end
    toc
    save('output','params','Likelihood')

end
