function dY = HATmodel(t,Y,eta,BV,muV0,muV1,sigmaV,aH,aL, ...
               betaVH,betaVL,tauV,muH,betaH,tauH,gammaH1,gammaH2,betaL,tauL,gammaL,deltaL, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);


%% Compartments
VP = Y(1); VS = Y(2); VE = Y(3); VI = Y(4); VR = Y(5);
HS = Y(6); HE = Y(7);  HI1 = Y(8); HI2 = Y(9); HR = Y(10);
LS = Y(11); LE = Y(12); LI = Y(13); LR = Y(14);
HC = Y(15);

V = VS+VE+VI+VR;
H = HS+HE+HI1+HI2+HR;
L = LS+LE+LI+LR;
% $$$ % If LiveStock is included
% $$$ if LivStock == 1
% $$$     LS = Y(12); LE = Y(13); LI = Y(14); LR = Y(15);
% $$$     L = LS+LE+LI+LR;
% $$$     lambdaVL = betaVL*LI/L;
% $$$ else
% $$$     aL = 0;
% $$$     lambdaVL = 0;
% $$$ end

% If vector-control being implemented
if rho == 0
    mr = 0;
else
    %    aH = aH*(1-mortality_rate(rho,t,l,m));
    %    aL = aL*(1-mortality_rate(rho,t,l,m));
    mr = mortality_rate(rho,t,l,m);
end



phi1 = P1*P1PD*P1TP;
phi2 = P2*P2PD*P2TP;
BH = muH*H + ((1-phi2)*gammaH2 + phi2*(1-eps2)*p2*zeta2)*HI2;
muV = muV0*(1+muV1*V);

lambdaVH = betaVH*HI1/H;
lambdaVL = betaVL*LI/L;

%% Model equations:

% Tsetsi Equations
dVP = BV*V - eta*VP;
dVS = eta*V - (aH+aL+ sigmaV + muV+mr)*VS;
dVE = (aH*lambdaVH+aL*lambdaVL)*VS - (tauV+muV+mr)*VE;
dVI = tauV*VE - (muV+mr)*VI;
dVR = (aH*(1-lambdaVH)+aL*(1-lambdaVL)+sigmaV)*VS -(muV+mr)*VR;

% Human Equations
dHS = BH + deltaH*HR - aH*betaVH*betaH*VI*HS/H - muH*HS;
dHE = aH*betaVH*betaH*VI*HS/H - (tauH+muH)*HE;
dHI1 = tauH*HE  -(phi1*eps1*zeta1 +(1-phi1)*gammaH1 + muH)*HI1;
dHI2 = (1-phi1)*gammaH1*HI1 - (phi2*eps2*zeta2  + (1-phi2)*gammaH2 + phi2*(1-eps2)*p2*zeta2 + muH)*HI2;
dHR = phi1*eps1*zeta1*HI1 + phi2*eps2*zeta2*HI2 -(deltaH+muH)*HR;

dHC = tauH*HE + (1-phi1)*gammaH1*HI1 + phi1*eps1*zeta1*HI1 + phi2*eps2*zeta2*HI2;

% Animal Equations
dLS = deltaL*LR - aL*betaVL*betaL*VI*LS/L ;
dLE = aL*betaVL*betaL*VI*LS/L - tauL*LE;
dLI = tauL*LE - gammaL*LI;
dLR = gammaL*LI-deltaL*LR;

dY = vertcat(dVP,dVS,dVE,dVI,dVR,dHS,dHE,dHI1,dHI2,dHR,dLS,dLE,dLI,dLR,dHC);

end
