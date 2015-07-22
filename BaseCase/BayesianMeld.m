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
    [a,b,c] = runHATmodel(Par(nonzeroind(i),:));
    A(i) = a;
    B(i) = b;
    C(i) = c;
end

X = betarnd(5,4307,10000,1);
Y = betarnd(7,4307,10000,1);

ci1 = quantile(X,[0.025,0.975])
ci2 = quantile(Y,[0.025,0.975])

t = 1:SampleSize;
t0 = [1,SampleSize]
figure;
subplot(3,1,1)
plot(t,A,'ro',t0,[ci1(1),ci1(1)],'b',t0,[ci1(2),ci1(2)],'b')
subplot(3,1,2)
plot(t,B,'ro',t0,[ci2(1),ci2(1)],'b',t0,[ci2(2),ci2(2)],'b')
subplot(3,1,3)
plot(t,C,'ro',t0,[0.0105,0.0105],'b')

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
