function dY = HATmodel(t,Y,eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);

%HATmodel(t, Y, eta, BV, muV0, muV1, sigmaV, aH, betaVH, tauV, muH, betaH, tauH, gammaH1, gammaH2, P1, P1PD, P1TP, P2, P2PD, P2TP, eps1, eps2, p2, deltaH, zeta1, zeta2, rho, l, m)
%% Compartments
VP = Y(1); VS = Y(2); VE = Y(3); VI = Y(4); VR = Y(5);
HS = Y(6); HE = Y(7);  HI1 = Y(8); HI2 = Y(9); HR = Y(10); HC = Y(11);

V = VS+VE+VI+VR;
H = HS+HE+HI1+HI2+HR;


% For vector-control
if rho == 0
    mr = 0;
else
    mr = mortality_rate(rho,t,l,m);
end

% coverage for stage 1 and stage 11
phi1 = P1*P1PD*P1TP;
phi2 = P2*P2PD*P2TP;

% Birth-rate for humans and death rate of tsetse
BH = muH*H + ((1-phi2)*gammaH2 + phi2*(1-eps2)*p2*zeta2)*HI2;
muV = muV0*(1+muV1*V);

lambdaVH = betaVH*HI1/H;

%% Model equations:

% Tsetse Equations
dVP = BV*V - eta*VP;
dVS = eta*V - (aH+ sigmaV + muV+mr)*VS;
dVE = (aH*lambdaVH)*VS - (tauV+muV+mr)*VE;
dVI = tauV*VE - (muV+mr)*VI;
dVR = (aH*(1-lambdaVH)+sigmaV)*VS -(muV+mr)*VR;

% Human Equations
dHS = BH + deltaH*HR - aH*betaVH*betaH*VI*HS/H - muH*HS;
dHE = aH*betaVH*betaH*VI*HS/H - (tauH+muH)*HE;
dHI1 = tauH*HE  -(phi1*eps1*zeta1 +(1-phi1)*gammaH1 + muH)*HI1;
dHI2 = (1-phi1)*gammaH1*HI1 - (phi2*eps2*zeta2  + (1-phi2)*gammaH2 + phi2*(1-eps2)*p2*zeta2 + muH)*HI2;
dHR = phi1*eps1*zeta1*HI1 + phi2*eps2*zeta2*HI2 -(deltaH+muH)*HR;
dHC = tauH*HE+(1-phi1)*gammaH1*HI1+phi1*eps1*zeta1*HI1 + phi2*eps2*zeta2*HI2;

dY = vertcat(dVP,dVS,dVE,dVI,dVR,dHS,dHE,dHI1,dHI2,dHR,dHC);


end
