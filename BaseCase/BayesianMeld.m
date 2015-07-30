% Par = [];
% Lik = [];
% N = 9  % No of independent sample runs that I ran
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

Data =    [3.0000    2.0000    9.5300
           4.0000    8.0000         0
           7.0000   13.0000         0
           3.0000    4.0000         0];

SampSize = [1488        1488        1634
            4514        4514           0
            7708        7708           0
            7788        7788           0];




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

% plot prevalences
t = 1:4;
fig1 = figure;
subplot(1,2,1)
plot(t,A(:,1:4),'o');
title('Stage I')
ax = gca;
ax.XTick = [1,2,3,4]
ax.XTickLabel = {'2008', '2010','2012',' 2013'}
subplot(1,2,2)
plot(t,A(:,5:8),'o');
title('Stage II')
ax = gca;
ax.XTick = [1,2,3,4]
ax.XTickLabel = {'2008', '2010','2012',' 2013'}

fig2 = figure;
plot(A(:,9),'o')
title('Vector Prevalence (2008)')


% Posterior
total = 500;

j = 1;
while j < total
    k =randi(length(weights),1);
    if(rand% plot prevalences
t = 1:4;
fig1 = figure;
subplot(1,2,1)
plot(t,A(:,1:4));
hold on;
plot(t,Data(:,1)./SampSize(:,1),'o','linewidth',2)
hold off;
title('Stage I')
ax = gca;
ax.XTick = [1,2,3,4]
ax.XTickLabel = {'2008', '2010','2012',' 2013'}
subplot(1,2,2)
plot(t,A(:,5:8));
hold on;
plot(t,Data(:,1)./SampSize(:,1),'o','linewidth',2)
hold off;
title('Stage II')
ax = gca;
ax.XTick = [1,2,3,4]
ax.XTickLabel = {'2008', '2010','2012',' 2013'}

fig2 = figure;
plot(A(:,9))
title('Vector Prevalence (2008)')


% Posterior
total = 500;

j = 1;
while j < total
    k =randi(length(weights),1);
    if(rand <weights(k))
    param = Par(k,:);
    posterior(j,:)=param;
     j = j+1;
    end
end


% Plot posterior

fig3 = figure;
for i = 1:4
    subplot(2,2,i)
    hist(posterior(:,i),20)
end
