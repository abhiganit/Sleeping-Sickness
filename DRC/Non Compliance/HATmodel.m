function dY = HATmodel(t,Y,eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,k,nuH,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);


%% Compartments
VS = Y(1); VE = Y(2); VI = Y(3); VR = Y(4);
HS = Y(5); HE = Y(6); HI1 = Y(7);  HA = Y(8); HI2 = Y(9); HR = Y(10); HC = Y(11);

V = VS+VE+VI+VR;
H = HS+HE+HI1+HA+HI2+HR;


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
muV = BV;

lambdaVH = betaVH*(HI1+k*HA)/H;

%% Model equations:

% Tsetsi Equations
dVS = BV - (aH+ sigmaV + muV+mr)*VS;
dVE = (aH*lambdaVH)*VS - (tauV+muV+mr)*VE;
dVI = tauV*VE - (muV+mr)*VI;
dVR = (aH*(1-lambdaVH)+sigmaV)*VS -(muV+mr)*VR;

% Human Equations
dHS = BH + deltaH*HR - aH*betaH*VI*HS - muH*HS;
dHE = aH*betaH*VI*HS - (tauH+muH)*HE;
dHI1 = nuH*tauH*HE  -(gammaH1 + muH)*HI1;
dHA  = (1-nuH)*tauH*HE - gammaH1*HA -muH*HA;
dHI2 = gammaH1*(HI1+HA) - (phi2*eps2*zeta2  + (1-phi2)*gammaH2 + phi2*(1-eps2)*p2*zeta2 + muH)*HI2;
dHR =  phi2*eps2*zeta2*HI2 + (gammaH1*gammaH2/(gammaH1+gammaH2))*HA-(deltaH+muH)*HR;
dHC = phi2*eps2*zeta2*HI2;

dY = vertcat(dVS,dVE,dVI,dVR,dHS,dHE,dHI1,dHA,dHI2,dHR,dHC);


end
