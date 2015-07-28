
function[] = modelfit()
%% Sampling.

    tic;

    X = betarnd(3,1488,10000,1);
    Y = betarnd(2,1488,10000,1);

    ci1 = quantile(X,[0.025,0.975])
    ci2 = quantile(Y,[0.025,0.975])


    N = 50000;
    % N samples from priors
    params = zeros(N,5);
    tic
    parfor i = 1:N
        betaVH = 0.1+0.5*rand;
        betaH = rand;
        zeta = 1.37*rand;
        P1 = rand;
        rho = 365*0.15*rand;

        params(i,:)  = [betaVH,betaH,zeta,P1,rho];


        [S1,S2,T,ST1,ST2,SV01,SV02,SV11,SV12] = ...
            runHATmodel(params(i,:));
        %        [S1,S2,T] = runHATmodel(params(i,:));

        A = betapdf(S1,3,1488);
        B = betapdf(S2,2,1488);
        C = betapdf(T,9.53,1634);
        %        C = betapdf(T,460,50000);
        D = betapdf(ST1,4,4514);
        E = betapdf(ST2,8,4514);
        F = betapdf(SV01,7,7708);
        G = betapdf(SV02,13,7708);
        H = betapdf(SV11,3,7788);
        I = betapdf(SV12,4,7788);

        if (S1<=ci1(1)) || (S1>=ci1(2)) || (S2<=ci2(1)) || (S2>=ci2(2)) || T > 0.010
	  Likelihood(i) = 0;
        else
	  Likelihood(i) = A*B*C*D*E*F*G*H*I;
        end


    end
    toc
    save('output','params','Likelihood')

end
