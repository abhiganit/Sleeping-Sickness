function[HI]  = runHATintervention()

%% Notes
%  1. BASE CASE: LivStock and AsymCarrier variable  should be set to
%     zero.
%  2. ANIMAL RESERVOIR: To explore role of animal reservoir, make
%     LivStock variable = 1 and AsymCarrier = 0.
%  3. ASYMTOMATIC CARRIERS: For asymptomatic carriers, set
%     AsymCarrier = 1 and LivStock = 0.

%% Interventions implemented



tic;


%% Tsetse Parameters (All rates are in years)
eta = 365./20;                       % 1/eta: pupae stage duration
BV = 365*0.05;                       % Tsetse constant birth rate
muV0 = 365*0.025;                    % Tsetse death rate without competition
muV1 = 0.0002;                       % Death rate competition parameters
sigmaV = 365.;                       % 1/sigmaV: Susceptibility period in Tsetse
aH = 365*0.075;                      % Tsetse human biting rate
betaVH = 2.515;                      % Tran. prob. from humans to Tsetse
tauV = 365./25;                      % 1/tauV: incubation period in tsetse
V = 328;                            % Tsetse population size (carrying capacity)

%% Human Parameters (All rates are in years)
muH = 1/55;                          % Human natural death rate
betaH = 0.62;                        % Trans. prob. from tsetse to humans
tauH = 365./12;                      % 1/tauH:incubation period in humans
gammaH1 = 365*0.0028;                % 1/gammaH1: 1st stage infectious period in humans
gammaH2 = 365*0.0029;                % 1/gammaH2: 2nd stage infectious period in humans
H = 300;                             % Human population size

%% Human Treatment Parameters
P1 = 0.68;                            % Prob. a stage I individual gets CATT test
P1PD = 0.87;                         % Sensitivity of test for stage I patientP
P1TP = 1;                            % Prob. of treatment after testing +ve for stage I patient
P2 = 1;                              % Prob. a stage II individual gets CATT test
P2PD = 0.87;                         % Sensitivity of test for stage II patient
P2TP = 1;                            % Prob. of treatment after testing +ve for stage II patient
eps1 = 0.94;                         % Stage I treatment efficacy
eps2 = 0.965;                        % Stage II treatment efficacy
p2 = 0.007;                          % Prob of treatment failure mortality in stage II patient
deltaH = 365./50;                    % 1/deltaH: immune period in  humans after treatment
zeta1 = 36.5;
zeta2 = 36.5;


%% Vector Treatment Parameters
rho = 0.0;                           % constrant mortality rate
l = 0;                               % l months at highest capacity
m = 0;                               % next m months of linear decline



%% Solve Model

% Initial conditions
V0 = [0.01*V,0.99*V, 0, 0.01*V,0];   % (Vp,Vs,Ve,Vi,Vr)
H0 = [H,0,0,0,0];                  % (Hs,He,Hc,Hi1,Hi2,Hr)
y0 = horzcat(V0,H0);



% Time span
%tspan = linspace(0,1000,1000) ;        % In years (for asym. case 1000->5000)
tspan = [0,1000];

% ODE solver (run model to equilibrium)
[t,y] = ode45(@HATmodel,tspan, y0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);


%% Running model with Interventions
ye0 = y(end,:); % Initializing from equilibrium
intervention = 2;

if intervention == 1  % increased treatment coverage
    t = [];
    y = [];
    tint = linspace(0,1,360);
    for j = 1:20
        [tl,yl] = ode45(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
        t = vertcat(t,tl);
        y = vertcat(y,yl);
        tint = linspace(tint(end)+10/360,tint(end)+1,1000);
        ye0  = yl(end,:);
        ye0(10) = ye0(10)+0.8*ye0(8);
        ye0(8) = 0.2*ye0(8);

    end
elseif intervention == 2 %
    tint = linspace(0,5,5*360);
    rho = 365*0.01;
    l = 3;
    m = 3;
    [ti,yi] = ode45(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
    tint = linspace(5,20,15*360);
    ye0 = yi(end,:);
    rho = 0;
    [tpi,ypi] = ode45(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
               betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
               P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
   t = vertcat(ti,tpi);
   y = vertcat(yi,ypi)
end





%% Plots
Compartments = {'VP','VS','VE','VI','VR','HS','HE','HI1','HI2','HR'}

for i = 1:length(Compartments)
    eval([Compartments{i} '= y(:,i);']);
end

HI = HI1(end) + HI2(end);

close all
figure;
subplot(2,2,1)
plot(t,VP,'g',t,VS,'linewidth',2);
% xlim([0,0.02]);
title('V_P & V_S')
xlabel('Years')
ylabel('No. of tsetse')
legend('V_P','V_S')
legend('boxoff','southeast')
box('off')
subplot(2,2,2)
plot(t, VE,'m','linewidth',2);
% xlim([0,1]);
title('V_E')
xlabel('Years')
box('off')
subplot(2,2,3)
plot(t, VI,'r','linewidth',2);
% xlim([0,1]);
title('V_I')
xlabel('Years')
ylabel('No. of tsetse')
box('off')
subplot(2,2,4)
plot(t , VR,'c','linewidth',2);
% xlim([0,1]);
title('V_R')
xlabel('Years')
box('off')

figure;
subplot(2,2,1)
plot(t,HS, 'g','linewidth',2);
%xlim([0,80]);
xlabel('Years')
ylabel('No. of people')
title('H_S')
box('off')
subplot(2,2,2)
plot(t,HE, 'm','linewidth',2);
%xlim([0,80]);
title('H_E')
box('off')
subplot(2,2,3)
plot(t, HI1,'r',t,HI2, 'linewidth',2);
xlabel('Years')
ylabel('No. of people')
%xlim([0,80]);
title('H_{I_1} & H_{I_2}')
legend('H_{I_1}','H_{I_2}')
legend('boxoff','southeast')
box('off')
subplot(2,2,4)
plot(t,HR,'c','linewidth',2);
title('H_{R}')
%xlim([0,80]);
xlabel('Years')
box('off')

toc

end
