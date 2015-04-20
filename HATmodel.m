function dY = HATmodel(t,Y,eta,BV,muV0,muV1,sigmaV,aH, ...
               aL,betaVH,betaVL,tauV,muH,betaH,tauH,nuH,k1,k2,gammaHC,gammaH1,gammaH2, ...
               betaL,tauL,gammaL,deltaL,P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,LivStock,alpha)


%% Compartments
VP = Y(1); VS = Y(2); VE = Y(3); VI = Y(4); VR = Y(5);
HS = Y(6); HE = Y(7); HC = Y(8); HI1 = Y(9); HI2 = Y(10); HR = Y(11);

V = VS+VE+VI+VR;
H = HS+HE+HC+HI1+HI2+HR;

% If LiveStock is included
if LivStock == 1
    LS = Y(12); LE = Y(13); LI = Y(14); LR = Y(15);
    L = LS+LE+LI+LR;
    lambdaVL = betaVL*LI/L;
else
    aL = 0;
    lambdaVL = 0;
end

% Composite parameters
phi1 = P1*P1PD*P1TP;
phi2 = P2*P2PD*P2TP;
BH = muH*H + ((1-phi2)*gammaH2 + phi2*(1-eps2)*p2)*HI2;
muV = muV0*(1+muV1*V);
lambdaVH = betaVH*(HI1+k1*HC+k2*HI2)/H;

%% Model equations:

% Tsetsi Equations
dVP = BV*V - eta*VP;
dVS = eta*V - (aH+ aL+ sigmaV + muV)*VS;
dVE = (aH*lambdaVH +aL*lambdaVL)*VS - (tauV+muV)*VE;
dVI = tauV*VE - muV*VI;
dVR = (aH*(1-lambdaVH)+aL*(1-lambdaVL)+sigmaV)*VS -muV*VR;

% Human Equations
dHS = BH + deltaH*HR - aH*betaH*VI*HS/H - muH*HS;
dHE = aH*betaH*VI*HS/H - (tauH+muH)*HE;
dHC = ((1-nuH)*tauH*HE - (gammaHC+muH)*HC)*alpha;
dHI1 = nuH*tauH*HE  -(phi1*eps1 +(1-phi1)*gammaH1 + muH)*HI1;
dHI2 = (1-phi1)*gammaH1*HI1 - (phi2*eps2  + (1-phi2)*gammaH2 + phi2*(1-eps2)*p2 + muH)*HI2;
dHR = phi1*eps1*HI1 + phi2*eps2*HI2 + alpha*gammaHC*HC -(deltaH+muH)*HR;

% if LiveStock is included
if LivStock == 1
    dLS = deltaL*LR - aL*betaL*VI*LS/L ;
    dLE = aL*betaL*VI*LS/L - tauL*LE;
    dLI = tauL*LE - gammaL*LI;
    dLR = gammaL*LI-deltaL*LR;

    dY = vertcat(dVP,dVS,dVE,dVI,dVR,dHS,dHE,dHC,dHI1,dHI2,dHR,dLS,dLE,dLI,dLR);
else
    dY = vertcat(dVP,dVS,dVE,dVI,dVR,dHS,dHE,dHC,dHI1,dHI2,dHR);
end

end
