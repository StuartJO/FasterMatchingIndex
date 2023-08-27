a = zeros(100);
D = rand(100);
D = triu(D,1)+triu(D,1)';

t1 = zeros(1,20);
t2 = t1;
t = t2;

for i = 1:20

m = 100*i;
tic
[B,b] = gen_model_mult(a,{D},m,'matching',{'exponential','powerlaw'},-2,.4);
t1(i) = toc;

tic
[B_,b_] = gen_model_mult_old(a,{D},m,'matching',{'exponential','powerlaw'},-2,.4);
t2(i) = toc;

t(i) = t2/t1;

end

plot(t)

save('MatchingGen2.mat')