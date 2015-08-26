function[out]  = runHATmodel(x)

%tic;

%% Tsetse Parameters (All rates are in years)
eta = 365./20;                       % 1/eta: pupae stage duration
BV = 365*0.05;                       % Tsetse constant birth rate
muV0 = 365*0.030;                    % Tsetse death rate without competition
muV1 = 0.0002;                       % Death rate competition parameters
sigmaV = 365.;                       % 1/sigmaV: Susceptibility period in Tsetse
aH = 365/3;                      % Tsetse human biting rate
aL = 365/3;                      % Tsetse Livestock biting rate
betaVH = x(1);                       % Tran. prob. from humans to
betaVL = 1-x(1);                      % Tran. prob. from Livestock to Tsetse                                     % Tsetse
tauV = 365./25;                      % 1/tauV: incubation period in tsetse
V = 5000;                             % Tsetse population size (carrying capacity)

%% Human Parameters (All rates are in years)
muH = 1/59;                          % Human natural death rate
betaH = x(2);                        % Trans. prob. from tsetse to humans
tauH = 365./12;                      % 1/tauH:incubation period in humans
gammaH1 = 365/526; %*0.0019;                % 1/gammaH1: 1st stage infectious period in humans
gammaH2 = 365/252; %*0.0020;                % 1/gammaH2: 2nd stage infectious period in humans
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


%% Animal Parameters
betaL = x(5);                        % Trans. prob. from tsetse to livestock
tauL = 365./12;                      % 1/tauL: incubation period in livestock
gammaL = 365./50;                    % 1/gammaL: infectious period in livestock
deltaL = 365./50;                    % 1/deltaL: immune period in livestock
L = 50;                              % Livestock population size


%% Vector Treatment Parameters
rho = 0.0;                           % constrant mortality rate
l = 0;                               % l months at highest capacity
m = 0;                               % next m months of linear decline


%% Solve Model

% Initial conditions
V0 = [BV*V/eta,0.99*V, 0, 0.01*V,0];   % (Vp,Vs,Ve,Vi,Vr)
H0 = [H,0,0,0,0];                  % (Hs,He,Hc,Hi1,Hi2,Hr)
L0 = [L,0,0,0];

y0 = horzcat(V0,H0,L0);

% Time span
tspan = [0,1000];
% ODE solver (solve model to equilibrium)
[t,y] = ode23s(@HATmodel,tspan, y0, [], eta,BV,muV0,muV1,sigmaV,aH,aL, ...
               betaVH,betaVL,tauV,muH,betaH,tauH,gammaH1,gammaH2,betaL,tauL,gammaL,deltaL, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);


S1(1) = (y(end,8))/sum(y(end,6:10)); % 2008
S2(1) = (y(end,9))/sum(y(end,6:10)); % 2008
V1 = y(end,4)/sum(y(end,2:5));       % 2008
L1 = y(end,13)/sum(y(end,11:end));
% In year 2008, there was active surveillance in Boffa East
% mainland with  coverage= (attendance*sensitivity)
y01 = y(end,:);
cov = 0.1026*0.87;
y01(10) = y01(10) + cov*(y01(8)+y01(9));
y01(8) = (1-cov)*y01(8); y01(9) = (1-cov)*y01(9);

tspan1 = linspace(0,2,3);  % (2008-2010) Run from dec 2007 (0), dec (2008) (1),
                           % dec (2009) (2)
[t1,y1] = ode23s(@HATmodel,tspan1, y01, [], eta,BV,muV0,muV1,sigmaV,aH,aL, ...
               betaVH,betaVL,tauV,muH,betaH,tauH,gammaH1,gammaH2,betaL,tauL,gammaL,deltaL, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);



S1(2) =(y1(end,8))/sum(y1(end,6:10)); % 2010
S2(2) =( y1(end,9))/sum(y1(end,6:10)); % 2010

% (2010-2012)
y02 = y1(end,:);
cov = 0.3119*0.87;
y02(10) = y02(10) + cov*(y02(8)+y02(9));
y02(8) = (1-cov)*y02(8); y02(9) = (1-cov)*y02(9);

tspan2 = linspace(0,2,3); % (2010-2012): Run from dec (2009) (0) to dec (2010) (1) to dec
% (2011) (2)
[t2,y2] = ode23s(@HATmodel,tspan2, y02, [], eta,BV,muV0,muV1,sigmaV,aH,aL, ...
               betaVH,betaVL,tauV,muH,betaH,tauH,gammaH1,gammaH2,betaL,tauL,gammaL,deltaL, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);


S1(3) = (y2(end,8))/sum(y2(end,6:10)); % 2012
S2(3) = (y2(end,9))/sum(y2(end,6:10)); % 2012
VecP = sum(y2(end,2:5));

% year 2012
rho = x(4);
l = 3; m = 3;

% End of 2011 (beg. 2012)
y03 = y2(end,:);
cov = 0.534*0.87;
y03(10) = y03(10) + cov*(y03(8)+y03(9));
y03(8) = (1-cov)*y03(8); y03(9) = (1-cov)*y03(9);

tspan3 = linspace(0,1,2); % (2012): Run from dec (2011) (0) to dec (2012) (1)
[t3,y3] = ode23s(@HATmodel,tspan3, y03, [], eta,BV,muV0,muV1,sigmaV,aH,aL, ...
               betaVH,betaVL,tauV,muH,betaH,tauH,gammaH1,gammaH2,betaL,tauL,gammaL,deltaL, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);

S1(4) = (y3(end,8))/sum(y3(end,6:10)); % 2013
S2(4) = (y3(end,9))/sum(y3(end,6:10)); % 2013
VecP1 = sum(y3(end,2:5));
out{1} = horzcat(S1,S2,V1,L1);
out{2} = y3(end,:);
out{3} = VecP1/VecP;
%toc

end
