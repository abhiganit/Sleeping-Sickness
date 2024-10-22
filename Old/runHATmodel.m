function[HI]  = runHATmodel()

%% Notes
%  1. BASE CASE: LivStock and AsymCarrier variable  should be set to
%     zero.
%  2. ANIMAL RESERVOIR: To explore role of animal reservoir, make
%     LivStock variable = 1 and AsymCarrier = 0.
%  3. ASYMTOMATIC CARRIERS: For asymptomatic carriers, set
%     AsymCarrier = 1 and LivStock = 0.


%% Animal reservoir
LivStock = 0;                        % Whether to include Livestock dynamics or not

%% Asymtomatic carriers
AsymCarrier = 0;
if AsymCarrier == 1
    alpha = 1;                       % Switch variable for asym. carr.
    nuH  = 0.6301;                   % Proportion of individuals becoming infectious.
else
    alpha = 0;                       % Switch variable for asym. carr.
    nuH = 1;                         % Proportion of individuals  becoming infectious


end

%% Tsetse Parameters (All rates are in years)
eta = 365./20;                       % 1/eta: pupae stage duration
BV = 365*0.05;                       % Tsetse constant birth rate
muV0 = 365*0.025;                    % Tsetse death rate without competition
muV1 = 0.0002;                       % Death rate competition parameters
sigmaV = 365.;                       % 1/sigmaV: Susceptibility period in Tsetse
aH = 365*0.075;                      % Tsetse human biting rate
aL = 365*0.175;                      % Tsetse Livestock biting rate
betaVH = 2.515;                      % Tran. prob. from humans to Tsetse
betaVL = 0.065;                      % Tran. prob. from Livestock to Tsetse
tauV = 365./25;                      % 1/tauV: incubation period in tsetse
V = 5000;                            % Tsetse population size (carrying capacity)

%% Human Parameters (All rates are in years)
muH = 1/55;                          % Human natural death rate
betaH = 0.62;                        % Trans. prob. from tsetse to humans
tauH = 365./12;                      % 1/tauH:incubation period in humans
gammaHC = 1/7;                       % 1/gammaHC: infectious period in asymptomatic carriers
gammaH1 = 365*0.0028;                % 1/gammaH1: 1st stage infectious period in humans
gammaH2 = 365*0.0029;                % 1/gammaH2: 2nd stage infectious period in humans
H = 300;                             % Human population size
k1 = 1;                              % Carriers contribution to transmission
k2 = 0;                              % Stage 2 contribution to transmission

%% Livestock Parameters (All rates are in years)
betaL = 0.62;                        % Trans. prob. from tsetse to livestock
tauL = 365./12;                      % 1/tauL: incubation period in livestock
gammaL = 365./50;                    % 1/gammaL: infectious period in livestock
deltaL = 365./50;                    % 1/deltaL: immune period in livestock
L = 50;                              % Livestock population size

%% Human Treatment Parameters
P1 = 0.68;                           % Prob. a stage I individual gets CATT test
P1PD = 0.87;                         % Sensitivity of test for stage I patient
P1TP = 1;                            % Prob. of treatment after testing +ve for stage I patient
P2 = 0.68;                           % Prob. a stage II individual gets CATT test
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
H0 = [H,0,0,0,0,0];                  % (Hs,He,Hc,Hi1,Hi2,Hr)
L0 = [L,0,0,0];                      % (Vs,Ve,Vi,Vr)

if LivStock == 1
    y0 = horzcat(V0,H0,L0);
else
    y0 = horzcat(V0,H0);
end


% Time span
tspan = linspace(0,1000,1000) ;      % In years


% ODE solver
[t,y] = ode45(@HATmodel,tspan, y0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
               aL,betaVH,betaVL,tauV,muH,betaH,tauH,nuH,k1,k2,gammaHC,gammaH1,gammaH2, ...
               betaL,tauL,gammaL,deltaL,P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
               eps2,p2,deltaH,zeta1,zeta2,LivStock,alpha,rho,l,m);



%% Plots
Compartments = {'VP','VS','VE','VI','VR','HS','HE','HC','HI1','HI2','HR','LS','LE','LI','LR'}

for i = 1:length(Compartments)-(1-LivStock)*4
    eval([Compartments{i} '= y(:,i);'])
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
plot(t,HC,'k',t, HI1,'r',t,HI2, 'linewidth',2);
xlabel('Years')
ylabel('No. of people')
%xlim([0,80]);
title('H_C, H_{I_1} & H_{I_2}')
legend('H_C','H_{I_1}','H_{I_2}')
legend('boxoff','southeast')
box('off')
subplot(2,2,4)
plot(t,HR,'c','linewidth',2);
title('H_{R}')
%xlim([0,80]);
xlabel('Years')
box('off')

if LivStock ==1
    figure;
    subplot(2,2,1)
    plot(t,LS,'g','linewidth',2);
    % xlim([0,0.2]);
    xlabel('Years')
    ylabel('No. of animals')
    title('L_S')
    box('off')
    subplot(2,2,2)
    plot(t, LE,'m','linewidth',2);
    % xlim([0,1]);
    xlabel('Years')
    title('L_E')
    box('off')
    subplot(2,2,3)
    plot(t, LI,'r','linewidth',2);
    % xlim([0,1]);
    xlabel('Years')
    ylabel('No. of animals')
    box('off')
    title('L_I')
    subplot(2,2,4)
    plot(t , LR,'c','linewidth',2);
    xlabel('Year')
    % xlim([0,1]);
    title('L_R')
    box('off')
end



end
