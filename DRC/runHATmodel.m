function[out]  = runHATmodel(x)

%tic;

%% Tsetse Parameters (All rates are in years)
eta = 365./20;                       % 1/eta: pupae stage duration
BV = 365*0.05;                       % Tsetse constant birth rate
muV0 = 365*0.030;                    % Tsetse death rate without competition
muV1 = 0.0002;                       % Death rate competition parameters
sigmaV = 365.;                       % 1/sigmaV: Susceptibility period in Tsetse
aH = 365/3;                          % Tsetse biting rate
betaVH = x(1);                       % Tran. prob. from humans to Tsetse
tauV = 365./25;                      % 1/tauV: incubation period in tsetse
V = 5000;                             % Tsetse population size (carrying capacity)
k = x(2);
%% Human Parameters (All rates are in years)
nuH = x(3);
muH = 1/59;                          % Human natural death rate
betaH = x(4);                        % Trans. prob. from tsetse to humans
tauH = 365./12;                      % 1/tauH:incubation period in humans
gammaH1 = 365/526;                   % 1/gammaH1: 1st stage infectious period in humans
gammaH2 = 365/252;                   % 1/gammaH2: 2nd stage infectious period in humans
H = 300;                             % Human population size

%% Human Treatment Parameters
P1 = 0.0;                            % Prob. a stage I individual gets CATT test
P1PD = 0.87;                         % Sensitivity of test for stage I patientP
P1TP = 1;                            % Prob. of treatment after testing +ve for stage I patient
P2 = 1;                              % Prob. a stage II individual gets CATT test
P2PD = 0.87;                         % Sensitivity of test for stage II patient
P2TP = 1;                            % Prob. of treatment after testing +ve for stage II patient
eps1 = 0.94;                         % Stage I treatment efficacy
eps2 = 0.965;                        % Stage II treatment efficacy
p2 = 0.007;                          % Prob of treatment failure mortality in stage II patient
deltaH = 365./50;                    % 1/deltaH: immune period in  humans after treatment
zeta1 = 1/1 ;
zeta2 = 1/x(5);

%% Vector Treatment Parameters
rho = 0.0;                           % constrant mortality rate
l = 0;                               % l months at highest capacity
m = 0;                               % next m months of linear decline


%% Solve Model
% Use equilibrium condition to set up initial conditions
D0 = 20;
phi2 = P2*P2PD*P2TP;
muV = muV0*(1+muV1*V);
mr =0;
HE  = D0/(nuH*tauH);
HI1 = nuH*tauH*HE/(gammaH1+muH);
HA  = (1-nuH)*tauH*HE/(gammaH1*gammaH2/(gammaH1+gammaH2)+muH);
HI2 = gammaH1*HI1/(phi2*eps2*zeta2  + (1-phi2)*gammaH2 + phi2*(1- ...
                                                  eps2)*p2*zeta2 + ...
                   muH);

HR  = phi2*eps2*zeta2*HI2/(deltaH+muH);
VS = BV*V/(aH+sigmaV+muV);
lambdaVH = betaVH*(HI1+k*HA)/H;
VE = aH*lambdaVH*VS/(tauV+muV);
VI = tauV*VE/muV;
VR = (aH*(1-lambdaVH)+sigmaV)*VS/(muV+mr);
HS  = (tauH+muH)*HE*H/(aH*betaH*VI);

% Initial conditions
V0 = [VS,VE,VI,VR];
H0 = [HS,HE,HI1,HA,HI2,HR,0];
% V0 = [BV*V/eta,0.99*V, 0, 0.01*V,0];   % (Vp,Vs,Ve,Vi,Vr)
% H0 = [H,0,0,0,0,0];                  % (Hs,He,Hi1,Hi2,Hr,Hc)
y0 = horzcat(V0,H0);

% Time span
tspan = [0,100];
% ODE solver (solve model to equilibrium)
[t,y] = ode23s(@HATmodel,tspan, y0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,k,nuH,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m)

out = y(end,end);

end
