
load('random200_data4topomdl.mat')

Nnodes = 100;
seed = zeros(Nnodes);

% Make a pretend distance matrix. We don't need a real one because this is
% a test of speed. 
D = A_dist;

% Number of edges to generate. I chose 650 because that is a similar number
% to real world brain networks using a 100 node parcellation in one
% hemisphere

Nedges = nnz(adjs{1})/2;

% Number of iterations to perform
iters = 100;
time_new = zeros(iters,1);
time_old = zeros(iters,1);

eta = -0.0663620917068697;
gam = 0.268238489046537;

A = adjs{1};
A_vals{1} = sum(A,2);
A_vals{2} = clustering_coef_bu(A);
A_vals{3} = betweenness_bin(A)';
A_vals{4} = D(triu(A,1) > 0);

maxKS_new = zeros(iters,1);
maxKS_old = zeros(iters,1);

% n = 50;
% seed_ = zeros(n);
% r = rand(n);
% R = triu(r,1)+triu(r,1)';
% Nedges_ = 1000;

%[B,b] = gen_model_add_old_(seed_,{R},Nedges_,'matching',{'exponential','powerlaw'},[-.21,1],[.21;3.96],'max');

%[B,b] = gen_model_add(seed_,{R},Nedges_,'matching',{'exponential','powerlaw'},[-.21,1],[.21;3.96],'max');

% [B,b] = gen_model_add_old_(seed,{D},Nedges,'matching',{'exponential','powerlaw'},[-.21,1],[.21;3.96],'max');
% 
% [B,b] = gen_model_add(seed,{D},Nedges,'matching',{'exponential','powerlaw'},[-.21,1],[.21;3.96],'max');
% 
parfor i = 1:iters
tic
B_ = gen_model_add(seed,{D},Nedges,'neighbors',{'exponential','powerlaw'},[-.21,1],[.21;3.96],'max');
time_new(i) = toc;
[maxKS_new(i)] = calc_maxKS(A_vals,D,B_); 
end

parfor i = 1:iters
tic
B = gen_model_add_old_(seed,{D},Nedges,'neighbors',{'exponential','powerlaw'},[-.21,1],[.21;3.96],'max');
time_old(i) = toc;    
[maxKS_old(i)] = calc_maxKS(A_vals,D,B);
end

figure
subplot(1,2,1)
boxplot([time_new time_old])
AverageSpeedUp = mean(time_old)/mean(time_new);
title(['Mean speed up = ',num2str(AverageSpeedUp)])
ylabel('Time in seconds')
xticklabels({'New code','Old code'})

subplot(1,2,2)
boxplot([maxKS_new maxKS_old])
xticklabels({'New code','Old code'})
ylabel('Model fit')