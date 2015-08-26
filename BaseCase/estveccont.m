function[output] = estveccont(x)


out = runHATmodel([0.1785 0.1547 0.6960  x]);
Data = [3,2,9.53;
        4,8, 0  ;
        7,13,0  ;
        3,4, 0 ];
SampSize = [1488,1488,1634;
            4514,4514,0;
            7708,7708,0;
            7788,7788,0];

% A = out{1};
% x = (A(1,4)-(Data(4,1)/SampSize(4,1))).^2;
% y = (A(1,8)-(Data(4,2)/SampSize(4,2))).^2;
%output = x+y;

output = abs(out{3}-0.5);
