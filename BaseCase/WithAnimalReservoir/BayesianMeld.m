%% Load parameter samples
Par = [];
Lik = [];
N = 3  % No of independent sample runs that I ran
for i = 1:N
    filename = sprintf('output%d.mat',i);
    load(filename)
    Par = vertcat(Par,params);
    Lik = horzcat(Lik,Likelihood);
end

%load output1
Lik = Likelihood;
Par = params;
M = length(Lik) % Total sample size
Data = [3,2,9.53;
        4,8, 0  ;
        7,13,0  ;
        3,4, 0 ];
SampSize = [1488,1488,1634;
            4514,4514,0;
            7708,7708,0;
            7788,7788,0];

%% Calculate 95% CI for data points
bnds01 = []; bnds02 = [];
for i = 1:4
    bnds01 = vertcat(bnds01,quantile(betarnd(Data(i,1),SampSize(i,1),10000,1), ...
                    [0.025,0.975]));
    bnds02 = vertcat(bnds02,quantile(betarnd(Data(i,2),SampSize(i,2),10000,1), ...
                    [0.025,0.975]));
end

bnds1 = abs(repmat(Data(:,1)./SampSize(:,1),1,2)-bnds01)
bnds2 = abs(repmat(Data(:,2)./SampSize(:,2),1,2)-bnds02)


%% Calculate Weights and discard zero Likelihood parameters
% parfor i = 1:M
%     weights(i) = Lik(i)/sum(Lik);
% end

nonzeroind  = find(Lik ~=0);
SampleSize = length(nonzeroind)

%%% Run the model again for non-zero samples
parfor i = 1:length(nonzeroind)
    out = runHATmodel(Par(nonzeroind(i),:));
    A(i,:) = out{1};
    P(i,:) = out{2};
end

%% Estimate the best fit for 2008 data
x = (A(:,1)-(Data(1,1)/SampSize(1,1))).^2;
y = (A(:,5)-(Data(1,2)/SampSize(1,2))).^2;
[a,b] = min((x+y))

bestpar = Par(nonzeroind(b),:)

%%% Next Estimate the vector control
% Run the model with Par(nonzeroind(b),:)) and estimate data from
% 2013
[p,fval] = fminsearch(@estveccont,0.04)
% p = 22.9924, fval = 2.0834e-06


%%% Run all of them again with new vector control
parfor i = 1:length(nonzeroind)
    out = runHATmodel([Par(nonzeroind(i),1:3),22.9924,Par(nonzeroind(i),5)]);
    A(i,:) = out{1};
    P(i,:) = out{2};
end

Lik = Lik(nonzeroind);
%weights = weights(nonzeroind);
Par  = Par(nonzeroind,:);
length(Par)


%% Plots

% Find subsample of Par that generates prevalence data that falls within
% confidence interval for year 2012 as well.
j = 1;
for i = 1:length(A)
      if A(i,4) >= bnds01(4,1) & A(i,4) <= bnds01(4,2) & A(i,8) >= bnds02(4,1)  & A(i,8) <= bnds02(4,2)
        B(j,:) = A(i,:); % Save respective prevaences
        out = A(i,:);
        Q(j,:) = P(i,:); % Save respective end points (initial
                         % conditions for predictions)
        X(j,:) = Par(i,:); % Saving respective parameters (for predictions)
        L(j) = betapdf(out(1),Data(1,1),SampSize(1,1))* ...
               betapdf(out(4),Data(4,1),SampSize(4,1))*...
               betapdf(out(5),Data(1,2),SampSize(1,2))* ...
               betapdf(out(8),Data(4,2),SampSize(4,2));% ...
                                                       % betapdf(out(9),Data(1,3),SampSize(1,3));
        % recalculating likelihood based on only 2 years
        j = j+1;
     end
end
length(B)
%save('initconds','B','Q','X');
%[a,b] = max(L);
out = runHATmodel([bestpar(1:3),22.9924,bestpar(5)]);
Best = out{1};



x = 1:4
[a,b] = min(B(:,1))
y1 = B(b,1:4)
[a,b] = max(B(:,1))
y2 = B(b,1:4)
X = [x,fliplr(x)]
Y = [y1,fliplr(y2)]

x = 1:4
[a,b] = min(B(:,5))
y1 = B(b,5:8)
[a,b] = max(B(:,5))
y2 = B(b,5:8)
X1 = [x,fliplr(x)]
Y1 = [y1,fliplr(y2)]


t = 1:4;
fig1 = figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
%plot(t,B(:,1:4),'Color',[0.75,0.75,0.75]);
fill(X,Y,[0.75,0.75,0.75],'EdgeColor','none')
box('off')
hold on;
plot(t,Best(1:4),'r','linewidth',1);
errorbar(t,Data(:,1)./SampSize(:,1), bnds1(:,1),bnds1(:,2),'ko', ...
         'linewidth',1.5)
box('off')
ylabel('Prevalence','FontSize',16)
xlabel('Years','FontSize',16)
hold off;
title('Stage I','FontSize',16)
ax = gca;
ax.XTick = [1,2,3,4]
ax.XTickLabel = {'2008', '2010','2012',' 2013'}
subplot(1,2,2)
%plot(t,B(:,5:8),'Color',[0.75,0.75,0.75]);
fill(X1,Y1,[0.75,0.75,0.75],'EdgeColor','none')
hold on;
plot(t,Best(5:8),'r','linewidth',1);
errorbar(t,Data(:,2)./SampSize(:,2), bnds2(:,1),bnds2(:,2),'ko', ...
         'linewidth',1.5)
xlabel('Years','FontSize',16)
hold off;
title('Stage II','FontSize',16)
ax = gca;
ax.XTick = [1,2,3,4]
ax.XTickLabel = {'2008', '2010','2012',' 2013'}

fig2 = figure;
plot(B(:,9),'o')
title('Vector Prevalence (2008)')

fig3 = figure;
plot(B(:,10),'o')
title('Vector Prevalence (2008)')







%% Re-sampling based on weights (meliding)
parfor i = 1:length(B)
    weights(i) = L(i)/sum(L);
end

total = 500;

j = 1;
while j < total+1
    k =randi(length(B),1);
    if(rand <weights(k))
    posterior(j,:)=X(k,:);
     j = j+1;
    end
end

% % % Posterior fits
% % parfor i = 1:total
% %     out = runHATmodel(posterior(i,:));
% %     A(i,:) = out;
% % end


% % Plot posterior

fig3 = figure;
for i = 1:5
    subplot(5,1,i)
    hist(posterior(:,i),20)
end
