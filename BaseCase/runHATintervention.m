function[t,y]  = runHATintervention(x,init,intervention)


%% Interventions implemented

tic;

%% Tsetse Parameters (All rates are in years)
eta = 365./20;                       % 1/eta: pupae stage duration
BV = 365*0.05;                       % Tsetse constant birth rate
muV0 = 365*0.025;                    % Tsetse death rate without competition
muV1 = 0.0002;                       % Death rate competition parameters
sigmaV = 365.;                       % 1/sigmaV: Susceptibility period in Tsetse
aH = 365*0.075;                      % Tsetse human biting rate
betaVH = x(1);                       % Tran. prob. from humans to Tsetse
tauV = 365./25;                      % 1/tauV: incubation period in tsetse
V = 1634;                             % Tsetse population size (carrying capacity)

%% Human Parameters (All rates are in years)
muH = 1/55;                          % Human natural death rate
betaH = x(2);                        % Trans. prob. from tsetse to humans
tauH = 365./12;                      % 1/tauH:incubation period in humans
gammaH1 = 365*0.0019;                % 1/gammaH1: 1st stage infectious period in humans
gammaH2 = 365*0.0020;                % 1/gammaH2: 2nd stage infectious period in humans
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
zeta2 = 1/x(3);



%% Solve Model forward
% Initial conditions
y0 = init;

% Coverage in year 2013

% active case-finding and treatment
cov = 0.5808*0.87;
y0(10) = y0(10) + cov*(y0(8)+y0(9));
y0(8) = (1-cov)*y0(8); y0(9) = (1-cov)*y0(9);
% vector-control
rho = x(4); % 0.01*365;
l = 3; m = 3;
% % run model for a year
tspan = linspace(0,1,360); % (2013): Run from dec (2012) (0) to dec (2013) (1)
[tp,yp] = ode23s(@HATmodel,tspan, y0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
                 betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
                 P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
                 eps2,p2,deltaH,zeta1,zeta2,rho,l,m);



%% Running model with Interventions
ye0 = yp(end,:); % Initializing from equilibrium
t = tp;
y = yp;
switch(intervention)
  case 'continue'
    rho = 0;
    tint = linspace(1,18,17*360);
    [tc,yc] = ode23s(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
                betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
                P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
                eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
    t = vertcat(t,tc);
    y = vertcat(y,yc);
  case 'yearly vector control'
    rho = x(4);
    tint = linspace(1,8,7*360);
    [tc,yc] = ode23s(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
                betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
                P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
                eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
    t = vertcat(t,tc);
    y = vertcat(y,yc);
    rho = 0;
    tint = linspace(8,20,12*360);
    [tc,yc] = ode23s(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
                betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
                P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
                eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
    t = vertcat(t,tc);
    y = vertcat(y,yc);
  otherwise
end


toc

end
