function[] = modelfit5()

%% Sampling.

    tic;

    X = betarnd(5,4307,10000,1);
    Y = betarnd(7,4307,10000,1);

    ci1 = quantile(X,[0.025,0.975])
    ci2 = quantile(Y,[0.025,0.975])


    N = 10000;
    % N samples from priors
    params = zeros(N,3);
    tic
    parfor i = 1:N
        betaVH = 0.1+0.5*rand;
        betaH = rand;
        zeta = 1.37*rand;
        params(i,:)  = [betaVH,betaH,zeta];


        [a,b,c] = runHATmodel(params(i,:));

        A = betapdf(a,5,4307);
        B = betapdf(b,7,4307);
        C = betapdf(c,10,1634);

        if (A<=ci1(1)) && (A>=ci1(2)) || (B<=ci2(1)) && (B>=ci2(2)) || C > 0.01
			 Likelihood(i) = 0;
        else
	  Likelihood(i) = A*B*C;
        end


    end
    toc
    save('output5','params','Likelihood')

end
