% Par = [];
% Lik = [];
% N = 2  % No of independent sample runs that I ran
% for i = 1:N
%     filename = sprintf('Sample%d.mat',i);
%     load(filename)
%     Par = vertcat(Par,params);
%     Lik = horzcat(Lik,Likelihood);
% end

load output
Lik = Likelihood;
Par = params;

M = length(Lik) % Total sample size

Data = [7,13,9.53;3,4,0];
SampSize = [7708,7708,1634;7788,7788,0];


bnds1 = []; bnds2 = [];
for i = 1:2
    bnds1 = vertcat(bnds1,quantile(betarnd(Data(i,1),SampSize(i,1),10000,1), ...
                    [0.025,0.975]));
    bnds2 = vertcat(bnds2,quantile(betarnd(Data(i,2),SampSize(i,2),10000,1), ...
                    [0.025,0.975]));
end

bnds1 = abs(repmat(Data(:,1)./SampSize(:,1),1,2)-bnds1)
bnds2 = abs(repmat(Data(:,2)./SampSize(:,2),1,2)-bnds2)


parfor i = 1:M
    weights(i) = Lik(i)/sum(Lik);
end

nonzeroind  = find(Lik ~=0);
SampleSize = length(nonzeroind)
parfor i = 1:length(nonzeroind)
    out = runHATmodel(Par(nonzeroind(i),:));
    A(i,:) = out;
end

weights = weights(weights~=0);
Par  = Par(weights~=0,:);
length(Par)

% find maximum likelihood
[a,b] = max(weights);

% plot prevalences
t = 1:2;
fig1 = figure;
subplot(1,2,1)
plot(t,A(:,1:2));
hold on;
plot(t,A(b,1:2),'k','linewidth',2);
errorbar(t,Data(:,1)./SampSize(:,1), bnds1(:,1),bnds1(:,2),'ko','linewidth',2)
hold off;
title('Stage I ratio')
ax = gca;
ax.XTick = [1,2]
ax.XTickLabel = {'2012',' 2013'}
subplot(1,2,2)
plot(t,A(:,3:4));
hold on;
plot(t,A(b,3:4),'k','linewidth',2);
errorbar(t,Data(:,2)./SampSize(:,2), bnds2(:,1),bnds2(:,2),'ko','linewidth',2)
hold off;
title('HAT')
ax = gca;
ax.XTick = [1,2]
ax.XTickLabel = {'2012',' 2013'}

fig2 = figure;
plot(A(:,5))
title('Vector Prevalence (2008)')


% % Posterior
% total = 500;

% j = 1;
% while j < total
%     k =randi(length(weights),1);
%     if(rand <weights(k))
%     param = Par(k,:);
%     posterior(j,:)=param;
%      j = j+1;
%     end
% end


% % Plot posterior

% fig3 = figure;
% for i = 1:4
%     subplot(2,2,i)
%     hist(posterior(:,i),20)
% endg
