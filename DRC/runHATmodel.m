function[Loglik]  = runHATmodel(x,Data)

%tic;
%global Data
%% Tsetse Parameters (All rates are in years)
eta = 365./20;                       % 1/eta: pupae stage duration
BV = 365*0.05;                       % Tsetse constant birth rate
muV0 = 365*0.030;                    % Tsetse death rate without competition
muV1 = 0.0002;                       % Death rate competition parameters
sigmaV = 365.;                       % 1/sigmaV: Susceptibility period in Tsetse
aH = 365/3;                          % Tsetse biting rate
betaVH = x(1);                       % Tran. prob. from humans to Tsetse
tauV = 365./25;                      % 1/tauV: incubation period in tsetse
V = 1;                             % Tsetse population size (carrying capacity)
k = x(2);
%% Human Parameters (All rates are in years)
nuH = x(3);
muH = 1/59;                          % Human natural death rate
betaH = x(4);                        % Trans. prob. from tsetse to humans
tauH = 365./12;                      % 1/tauH:incubation period in humans
gammaH1 = 365/526;                   % 1/gammaH1: 1st stage infectious period in humans
gammaH2 = 365/252;                   % 1/gammaH2: 2nd stage infectious period in humans
H = Data(1,5);                             % Human population size

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


scal = x(6);
%% Solve Model
% Use equilibrium condition to set up initial conditions
D0 = Data(1,3);  % Initial point
phi1 = Data(1,2);  % Covarage of active detection for inital year
phi2 = P2*P2PD*P2TP;
muV = BV;
mr =0;


B1 = BV/(aH+sigmaV+muV);
B2 = aH*betaVH*B1/(H*(tauV+muV));
B3 = (tauV/muV)*B2;

C1 = nuH*tauH/(gammaH1+muH);
C2 = (1-nuH)*tauH/((gammaH1*gammaH2/(gammaH1+gammaH2))+muH);
C3 = gammaH1*C1/(phi2*eps2*zeta2  + (1-phi2)*gammaH2 + phi2*(1- ...
                                                  eps2)*p2*zeta2 + ...
              muH) ;
C4 = phi2*eps2*zeta2*C3/(deltaH+muH);
% C5 = (gammaH1*gammaH2/(gammaH1+gammaH2))/(deltaH+muH);
% C6 = betaVH/H;
% C7 = (tauH+muH)/(aH*betaH);


HS = (tauH+muH)/(aH*betaH*B3*(C1+k*C2));


HE = (H-HS)/(1+C1+C2+C3+C4);
HI1 = C1*HE;
HA = C2*HE;
HI2 = C3*HE;
HR = C4*HE;
VS = B1;
VE = B2*(C1+k*C2)*HE;
VI = B3*(C1+k*C2)*HE;
VR = 1-(VS+VE+VI);

if H <= HS || VE+VI > 0.01
    Loglik = -Inf;
    return 
end
% Initial conditions
V0 = [VS,VE,VI,VR];
H0 = [HS,HE,HI1,HA,HI2,HR,0];
% V0 = [BV*V/eta,0.99*V, 0, 0.01*V,0];   % (Vp,Vs,Ve,Vi,Vr)
% H0 = [H,0,0,0,0,0];                  % (Hs,He,Hi1,Hi2,Hr,Hc)
y0 = horzcat(V0,H0);



% Time span
tspan = [0,1];
%tspan = linspace(0,1,365);
% ODE solver (solve model to equilibrium)
scal = Data(1,2);
T = []; Y = [];
for i = 1:4
    [t,y] = ode45(@HATmodel,tspan, y0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,k,nuH,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
    tspan = [i+(1/25),i+1];
    y0 = y(end,:);
    y0(10) = y0(10) + scal*P1PD*(y0(7)+y0(9));
    y0(7) = y0(7) - scal*P1PD*y0(7);
    y0(9) = y0(9) - scal*P1PD*y0(9);
    y0(11) =0;
    T = vertcat(T,t);
    Y = vertcat(Y,y);
end
%plot(T,Y)


tspan = [0,1];

for i = 1:12
    [t,y] = ode45(@HATmodel,tspan, y0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,k,nuH,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
    tspan = [i+(1/25),i+1];
    out{1}(i) = y(end,end);
    y0 = y(end,:);
    y0(10) = y0(10) + Data(i,2)*P1PD*(y0(7)+y0(9));
    y0(7) = y0(7) - Data(i,2)*P1PD*y0(7);
    y0(9) = y0(9) - Data(i,2)*P1PD*y0(9);
    out{2}(i) = Data(i,2)*P1PD*(y0(7)+y0(9));
    y0(11) =0;
    x1 = Data(i,3);
    n1 = Data(i,2)*sum(H0);
    pr1 = out{2}(i)/Data(1,5);
    x2 = Data(i,4);
    n2 = sum(H0);
    pr2 = out{1}(i)/Data(1,5);

    out{3}(i) = log(binopdf(x1,ceil(n1),pr1))+log(binopdf(x2,ceil(n2),pr2));
    %out{3}(i) = log(poisspdf(x1, out{2}(i)))+log(poisspdf(x2, out{1}(i)));
    % out{3}(i) = (binopdf(x1,ceil(n1),pr1))*(binopdf(x2,ceil(n2),pr2));
  

end

Loglik =  sum(out{3});

end
