
Par = [];
Lik = [];
N = 7  % No of independent sample runs that I ran
for i = 1:N
    filename = sprintf('Sample%d.mat',i);
    load(filename)
    Par = vertcat(Par,params);
    Lik = horzcat(Lik,Likelihood);
end

M = length(Lik); % Total sample size


parfor i = 1:M
    weights(i) = Lik(i)/sum(Lik);
end


nonzeroind  = find(Lik ~=0);

nonzeroind = find(Likelihood ~=0)

Par = params;

parfor i = 1:length(nonzeroind)
    [a,b,c] = runHATmodel(Par(nonzeroind(i),:));
    A(i) = a;
    B(i) = b;
    C(i) = c;
end



X = betarnd(5,4307,10000,1);
Y = betarnd(7,4307,10000,1);
Z = betarnd(10,1634,10000,1);

ci1 = quantile(X,[0.025,0.975])
ci2 = quantile(Y,[0.025,0.975])
ci3 = quantile(Z,[0.025,0.975])


  % weight = weight(weight~=0);
    %params = params(weight~=0,:);

    % {nter1 = Desired number of sample for poserior
     %  k1 = 1 ;

    %}

   % total = 500;
    %j = 1;
   % while j < total
    % k =randi(length(weight),1);

    % if(rand <weight(k))
    % Par =params(k,:);
     % posterior(j,:)=Par;
     % j = j+1
    % end
