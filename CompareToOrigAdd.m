
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

eta = -0.131390256116879;
gam = 0.909685257015132;
alpha = 3.25337137754728;

A = adjs{1};
A_vals{1} = sum(A,2);
A_vals{2} = clustering_coef_bu(A);
A_vals{3} = betweenness_bin(A)';
A_vals{4} = D(triu(A,1) > 0);

maxKS_new = zeros(iters,1);
maxKS_old = zeros(iters,1);

B_ = gen_model_add(A,{D},Nedges+1,'matching',{'exponential','powerlaw'},[eta,1],[gam;alpha],'max');

B = gen_model_add_normmax(A,D,Nedges+1,'matching',{'exponential','powerlaw','powerlaw'},eta,gam,[1 alpha],0);


for i = 1:iters
tic
B_ = gen_model_add(seed,{D},Nedges,'matching',{'exponential','powerlaw'},[eta,1],[gam;alpha],'max');
time_new(i) = toc;
[maxKS_new(i)] = calc_maxKS(A_vals,D,B_); 
end


for i = 1:iters
tic
[B,b] = gen_model_add_normmax(zeros(100),D,Nedges,'matching',{'exponential','powerlaw','powerlaw'},eta,gam,[1 alpha],0);
time_old(i) = toc;    
[maxKS_old(i)] = calc_maxKS(A_vals,D,B);
end

clear B B_ i

%save('additiveSpeedTestData.mat')

jitterOffset = .5;

figure('Position',[232 262 1040  583])
subplot(1,2,1)
%boxchart([time_new time_old])
cmap = lines(2);
jitter = (rand(length(time_new),1)-.5)*jitterOffset;
x=ones(length(time_new),1);
scatter(x+jitter,time_new,'filled','MarkerFaceColor',cmap(1,:),'MarkerFaceAlpha',0.6);
hold on
jitter = (rand(length(time_old),1)-.5)*jitterOffset;
x=ones(length(time_old),1)*2;
scatter(x+jitter,time_old,'filled','MarkerFaceColor',cmap(2,:),'MarkerFaceAlpha',0.6);

plotjitterOffset = jitterOffset/2;

plot([1-plotjitterOffset 1+plotjitterOffset],[mean(time_new) mean(time_new)],'LineWidth',2,'Color','k')
plot([2-plotjitterOffset 2+plotjitterOffset],[mean(time_old) mean(time_old)],'LineWidth',2,'Color','k')

AverageSpeedUp = mean(time_old)/mean(time_new);
title(['Mean speed up = ',num2str(AverageSpeedUp)])
ylabel('Time in seconds')
xticks(1:2)
xticklabels({'New code','Old code'})
set(gca,'FontSize',18)
xlim([.5 2.5])

subplot(1,2,2)
%boxchart([maxKS_new maxKS_old])
cmap = lines(2);
jitter = (rand(length(maxKS_new),1)-.5)*jitterOffset;
x=ones(length(maxKS_new),1);
scatter(x+jitter,maxKS_new,'filled','MarkerFaceColor',cmap(1,:),'MarkerFaceAlpha',0.6);
hold on
jitter = (rand(length(maxKS_old),1)-.5)*jitterOffset;
x=ones(length(maxKS_old),1)*2;
scatter(x+jitter,maxKS_old,'filled','MarkerFaceColor',cmap(2,:),'MarkerFaceAlpha',0.6);

plotjitterOffset = jitterOffset/2;

plot([1-plotjitterOffset 1+plotjitterOffset],[mean(maxKS_new) mean(maxKS_new)],'LineWidth',2,'Color','k')
plot([2-plotjitterOffset 2+plotjitterOffset],[mean(maxKS_old) mean(maxKS_old)],'LineWidth',2,'Color','k')

xticks(1:2)
xticklabels({'New code','Old code'})
ylabel('Model fit')
set(gca,'FontSize',18)
xlim([.5 2.5])