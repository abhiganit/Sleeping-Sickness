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
    if(rand <weights(k))
    param = Par(k,:);
    posterior(j,:)=param;
     j = j+1;
    end
end


% Plot posterior

fig3 = figure;
for i = 1:3
    subplot(3,1,i)
    hist(posterior(:,i),20)
end
