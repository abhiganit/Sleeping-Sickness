function[out] = calcprior(x)
out = 1;
A = zeros(6); B = [5,0.01,0.5,5,1.44, 1];
for i = 1:length(x);
out = unifpdf(x(i),A(i),B(i))*out;
end


