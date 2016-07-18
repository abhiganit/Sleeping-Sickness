
Runforward = {'continue','yearly vector control','yearly vector control with active case-finding','yearly vector control with alternate active case-finding'}

for s = 2:4
    red = [0.25,0.5,0.75]
    j = 1
    for r = red
        load initconds
        X(:,4) = (1-r)*22.9924;
        intervention = Runforward{s}
        Y = [];
        for i = 1:length(B)
            [t,y] = runHATintervention(X(i,:),Q(i,:),intervention,r);
            Y{i} = y;
        end

        Fail = [];
        for i = 1:length(B)
            if size(Y{i},1)<length(Y{1})
                Fail = [Fail,i];
            end
        end
        Y(Fail) = [];
        X(Fail,:) = [];
        Q(Fail,:) = [];
        B(Fail,:) = [];


        % Calculate total prevalences at each time point till 2030;
        Shai = (cellfun(@(y) (y(:,8)+y(:,9))'./sum(y(:,6:end)'),Y,'UniformOutput',false));
        TotPrev = cell2mat(Shai');
        VP = cell2mat((cellfun(@(y) y(:,4)'./sum(y(:,2:5)'), Y,'UniformOutput',false))');


        for i = 1:size(TotPrev,2)
            Tot = TotPrev(:,i);
            E1(i) = length(Tot(1000*Tot<1)); % Eradication criteria
            E2(i) = length(Tot(10000*Tot<1)); % Elimination as public
                                              % health c   riteria
            E3(i) = length(Tot(100000*Tot<1));
            P1(i) = E1(i)/length(Tot);
            P2(i) = E2(i)/length(Tot);
            P3(i) = E3(i)/length(Tot);
        end
        [a,b] = min(abs(t-8));
        ProbMat1(j) = P2(b);
        ProbMat2(j) = P3(b);
        ProbCurve{j} = P2;
        Time{j} = t;
        filename = sprintf('Sensitivity%d.mat',s);
        save(filename,'ProbCurve','ProbMat1','ProbMat2','Time')
        j = j+1;
    end
end

k = 1
for s = 2:4
    filename = sprintf('Sensitivity%d.mat',s)
    load(filename)
    MatProb1(k,:) = ProbMat1;
    MatProb2(k,:) = ProbMat2;
    k = k+1;
end
