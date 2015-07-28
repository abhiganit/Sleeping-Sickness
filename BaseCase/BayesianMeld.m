Par = [];
Lik = [];
N = 9  % No of independent sample runs that I ran
for i = 1:N
    filename = sprintf('Sample%d.mat',i);
    load(filename)
    Par = vertcat(Par,params);
    Lik = horzcat(Lik,Likelihood);
end

M = length(Lik) % Total sample size


parfor i = 1:M
    weights(i) = Lik(i)/sum(Lik);
end


nonzeroind  = find(Lik ~=0);
SampleSize = length(nonzeroind)
parfor i = 1:length(nonzeroind)
    [a,b,c,d,e,f,g,h,k] = runHATmodel(Par(nonzeroind(i),:));
    A(i,:) = [a,b,c,d,e,f,g,h,k]
end

weights = weights(weights~=0);
Par  = Par(weights~=0,:);

length(Par)
    % {nter1 = Desired number of sample for poserior
     %  k1 = 1 ;

    %}
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

figure;
for i = 1:3
    subplot(3,1,i)
    hist(posterior(:,i),20)
end
