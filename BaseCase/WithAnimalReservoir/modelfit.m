function[] = modelfit()
%% Sampling.

    tic;

    % Data & Sample: Years(2008/10/12/13). Columns:  Stage 1 | Stage 2
    % Vector Prevalence (0.049/50 and 0.8/300)
    c = clock;
    Time = fix(c);
    Time(4:end)

    Data = [3,2,9.53;
            4,8, 0.049  ;
            7,13,0  ;
            3,4, 0 ];
    SampSize = [1488,1488,1634;
                4514,4514,50;
                7708,7708,0;
                7788,7788,0];

    X = betarnd(3,1488,10000,1);
    Y = betarnd(2,1488,10000,1);
    ci1 = quantile(X,[0.025,0.975]) % S1, 2008
    ci2 = quantile(Y,[0.025,0.975]) % S2, 2008


    N = 1000000;
    % N samples from priors
    params = zeros(N,5);
    tic
    parfor j = 1:N
        betaVH = 0.1+0.5*rand;
        betaH = rand;
        zeta = 0.7*rand;
        rho = rand;
        betaL = rand;
        params(j,:)  = [betaVH,betaH,zeta,rho,betaL];

        out = runHATmodel(params(j,:));

        if (out{1}(1)<=ci1(1)) || (out{1}(1)>=ci1(2)) || (out{1}(5)<=ci2(1)) ...
                    || (out{1}(5)>=ci2(2)) || out{1}(9) > 0.010 || out{1}(9) < 10^-5 ||...
                out{1}(10) > 0.010
            Likelihood(j) = 0;
        else
            Lik1 = 1; Lik2 = 1;
            for i = 1:4
                Lik1 = betapdf(out{1}(i),Data(i,1),SampSize(i,1));
                Lik2 = betapdf(out{1}(i+4),Data(i,2),SampSize(i,2));
            end
            Likelihood(j) = Lik1*Lik2*betapdf(out{1}(9),Data(1,3),SampSize(1,3))*betapdf(out{1}(10),Data(2,3),SampSize(2,3));
        end
    end
    toc
    save('Sample','params','Likelihood')

    c = clock;
    Time = fix(c);
    Time(4:end)


end