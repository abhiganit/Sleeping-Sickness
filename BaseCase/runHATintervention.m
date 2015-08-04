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
rho = x(4);
l = 3; m = 3;
% run model for a year
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
    tint = linspace(1,8,7*360);
    [tc,yc] = ode23s(@HATmodel,tint, y0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
                betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
                P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
                eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
    t = vertcat(t,tc);
    y = vertcat(y,yc);
  case 'yearly vector control'

  otherwise
end


% if intervention == 1  % increased treatment coverage
%     t = [];
%     y = [];
%     tint = linspace(0,1,360);
%     for j = 1:7
%         [tl,yl] = ode45(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
%                betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
%                P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
%                eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
%         t = vertcat(t,tl);
%         y = vertcat(y,yl);
%         tint = linspace(tint(end)+10/360,tint(end)+1,1000);
%         ye0  = yl(end,:);
%         ye0(10) = ye0(10)+0.8*ye0(8);
%         ye0(8) = 0.2*ye0(8);

%     end
% elseif intervention == 2 %
%     tint = linspace(0,5,5*360);
%     rho = 365*0.01;
%     l = 3;
%     m = 3;
%     [ti,yi] = ode45(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
%                betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
%                P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
%                eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
%     tint = linspace(5,20,15*360);
%     ye0 = yi(end,:);
%     rho = 0;
%     [tpi,ypi] = ode45(@HATmodel,tint, ye0, [], eta,BV,muV0,muV1,sigmaV,aH, ...
%                betaVH,tauV,muH,betaH,tauH,gammaH1,gammaH2, ...
%                P1,P1PD,P1TP,P2,P2PD,P2TP,eps1, ...
%                eps2,p2,deltaH,zeta1,zeta2,rho,l,m);
%    t = vertcat(ti,tpi);
%    y = vertcat(yi,ypi)
% end





% %% Plots
% Compartments = {'VP','VS','VE','VI','VR','HS','HE','HI1','HI2','HR'}

% for i = 1:length(Compartments)
%     eval([Compartments{i} '= y(:,i);']);
% end

% HI = HI1(end) + HI2(end);

% close all
% figure;
% subplot(2,2,1)
% plot(t,VP,'g',t,VS,'linewidth',2);
% % xlim([0,0.02]);
% title('V_P & V_S')
% xlabel('Years')
% ylabel('No. of tsetse')
% legend('V_P','V_S')
% legend('boxoff','southeast')
% box('off')
% subplot(2,2,2)
% plot(t, VE,'m','linewidth',2);
% % xlim([0,1]);
% title('V_E')
% xlabel('Years')
% box('off')
% subplot(2,2,3)
% plot(t, VI,'r','linewidth',2);
% % xlim([0,1]);
% title('V_I')
% xlabel('Years')
% ylabel('No. of tsetse')
% box('off')
% subplot(2,2,4)
% plot(t , VR,'c','linewidth',2);
% % xlim([0,1]);
% title('V_R')
% xlabel('Years')
% box('off')

% figure;
% subplot(2,2,1)
% plot(t,HS, 'g','linewidth',2);
% %xlim([0,80]);
% xlabel('Years')
% ylabel('No. of people')
% title('H_S')
% box('off')
% subplot(2,2,2)
% plot(t,HE, 'm','linewidth',2);
% %xlim([0,80]);
% title('H_E')
% box('off')
% subplot(2,2,3)
% plot(t, HI1,'r',t,HI2, 'linewidth',2);
% xlabel('Years')
% ylabel('No. of people')
% %xlim([0,80]);
% title('H_{I_1} & H_{I_2}')
% legend('H_{I_1}','H_{I_2}')
% legend('boxoff','southeast')
% box('off')
% subplot(2,2,4)
% plot(t,HR,'c','linewidth',2);
% title('H_{R}')
% %xlim([0,80]);
% xlabel('Years')
% box('off')

toc

end
