%% Compare the matching index with different numbers of nodes and a consistent density
Iters = 10;

NetSizes = 100:100:1000;

nNetSizes = length(NetSizes);

time_old = zeros(nNetSizes,Iters);
time_new = zeros(nNetSizes,Iters);

% Set the desired density

den = .2;

for i = 1:nNetSizes
    Nnodes = NetSizes(i);
    % Make a random network with a density of around .2
    
    % We care about testing speed, so no need to make multiple versions of
    % a network    
    r = rand(Nnodes);
    r_ = triu(r,1)+triu(r,1)';
    A = double(r_>(1-den));
    
    % In MATLAB it seems like the first time a script/function is called
    % upon, it gets loaded into memory and this can take some time. Usually
    % this is so small you would never noticed but I do notice it here when
    % comparing timings and it throws them out somewhat! So on the first
    % loop we just run it once to get it into memory
    
    if i == 1
        m1 = matching(A); 
        m2 = matching_ind_und(A);
    end
    
    for j = 1:Iters
              
        tic
        m2 = matching_ind_und(A);
        time_old(i,j) = toc;

        tic
        m1 = matching(A);
        time_new(i,j) = toc;     
    end
    
    disp(['Completed nodes = ',num2str(Nnodes),', old method = ',num2str(mean(time_old(i,:))),' seconds, new method = ',num2str(mean(time_new(i,:))),' seconds'])
end

clear i j A r r_ Nnodes m1 m2

save('matchingSpeedTestData1.mat')

%

figure('Position',[228 523 1844 563])
subplot(1,3,1)
plot(NetSizes,mean(time_new,2),'LineWidth',4)
hold on
plot(NetSizes,mean(time_old,2),'LineWidth',4)
legend('New code','Old code','Location','northwest')
ylabel('Time in seconds')
xlabel('Number of nodes in network')
set(gca, 'FontSize', 18)
xlim([min(NetSizes) max(NetSizes)])

subplot(1,3,2)
plot(NetSizes,mean(time_new,2),'LineWidth',4)
hold on
plot(NetSizes,mean(time_old,2),'LineWidth',4)
legend('New code','Old code','Location','northwest')
ylabel('Time in seconds (log scale)')
xlabel('Number of nodes in network')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 18)
xlim([min(NetSizes) max(NetSizes)])

subplot(1,3,3)
plot(NetSizes,mean(time_old,2)./mean(time_new,2),'LineWidth',4)
ylabel('Speed up factor of new code')
xlabel('Number of nodes in network')
set(gca, 'FontSize', 18)
xlim([min(NetSizes) max(NetSizes)])

print('./images/MatchingDemo1.svg','-dsvg')

clear all

%% Compare the matching generative network model with the old and new code

% Set up a network

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

for i = 1:iters
    tic
    B_ = gen_model_mult(seed,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
    time_new(i) = toc;
    [maxKS_new(i)] = calc_maxKS(A_vals,D,B_); 
end

for i = 1:iters
    tic
    B = gen_model_mult_old(seed,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
    time_old(i) = toc;    
    [maxKS_old(i)] = calc_maxKS(A_vals,D,B);
end

clear B B_ i adjs

save('matchingSpeedTestData2.mat')

%

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

print('./images/MatchingDemo2.svg','-dsvg')

clear all

%% Compare the matching generative network model with the old and new code

NetSizes = [100 250 500 1000 2000];

iters = 10;
nNetSizes = length(NetSizes);

time_new = zeros(nNetSizes,iters);
time_old = zeros(nNetSizes,iters);

timecourse_new = cell(nNetSizes,1);
timecourse_old = cell(nNetSizes,1);

eta = -2;
gam = .4;

a = zeros(10);
D = rand(10);
D = triu(D,1)+triu(D,1)';

b1 = matching_gen_model_mult(a,{D},10,{'exponential','powerlaw'},eta,gam);
b2 = matching_gen_model_mult_old(a,{D},10,{'exponential','powerlaw'},eta,gam);

for i = 1:nNetSizes
    
    Nnodes = NetSizes(i);
    a = zeros(Nnodes);
    D = rand(Nnodes);
    D = triu(D,1)+triu(D,1)';
    Nedges = 2500;  
    
    tnew_temp = zeros(iters,Nedges);
    told_temp = zeros(iters,Nedges);
    
        for k = 1:iters

        tStart = tic;
        [~,~,tnew_temp(k,:)] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
        time_new(i,k) = toc(tStart);

        tStart = tic;
        [~,~,told_temp(k,:)] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
        time_old(i,k) = toc(tStart);
        end  
        
        timecourse_new{i} = tnew_temp;
        timecourse_old{i} = told_temp;
        disp(['Completed nodes = ',num2str(Nnodes),', edges = ',num2str(Nedges),', old method = ',num2str(mean(time_old(i,:))),' seconds, new method = ',num2str(mean(time_new(i,:))),' seconds'])

end

clear tStart tnew_temp told_temp a D i k b1 b2

save('matchingSpeedTestData3.mat')

%

figure('Position',[326 218 1541 652])

subplot(1,2,1)

% lgd = cell(1,nNetSizes*2);
lines_cmap = lines(nNetSizes);

for i = 1:nNetSizes
   newtimes = mean(timecourse_new{i});
   oldtimes = mean(timecourse_old{i});
   plot(cumsum(newtimes),'-','Color',lines_cmap(i,:),'LineWidth',2)
%   lgd{sub2ind([2 nNetSizes],1,i)} = ['New code (' num2str(NetSizes(i)),')'];
   hold on
   plot(cumsum(oldtimes),'--','Color',lines_cmap(i,:),'LineWidth',2)
%   lgd{sub2ind([2 nNetSizes],2,i)} = ['Old code (' num2str(NetSizes(i)),')'];
end
xlabel('Iteration/Edges')
ylabel('Cumulative time in seconds')
%legend(gca,lgd,'NumColumns',5,'Location','northeast','FontSize',12)
set(gca,'FontSize',18)
set(gca, 'YScale', 'log')
xlimits = xlim;
ylimits = ylim;

p(1) = plot(-1:-1:-100,-1:-1:-100,'-','Color','k');
hold on
p(2) = plot(-1:-1:-100,-1:-1:-100,'--','Color','k');
for i = 1:nNetSizes
   s(i) = scatter(-1,-1,'filled','MarkerEdgeColor',lines_cmap(i,:),'MarkerFaceColor',lines_cmap(i,:));
   NetSizesLgd{i} = num2str(NetSizes(i));
end
xlim(xlimits)
ylim(ylimits)
lg1 = legend(gca,p,{'New code','Old code'},'NumColumns',2,'Location','south','FontSize',18);
ah1=axes('position',get(gca,'position'),'visible','off');

lg2 = legend(ah1,s,NetSizesLgd,'NumColumns',5,'Location','south','FontSize',18);
xlim(xlimits)
ylim(ylimits)

lg1.Title.String = 'Code version';
lg2.Title.String = 'Network size (nodes)';

lg1_orig_pos = lg1.Position;

lg1.Position(2) = lg1_orig_pos(2)+lg1_orig_pos(4);

subplot(1,2,2)
lgd2 = cell(1,nNetSizes);
for i = 1:nNetSizes
    newtimes = mean(timecourse_new{i});
    oldtimes = mean(timecourse_old{i});
   Improvement = cumsum(oldtimes)./cumsum(newtimes);
   plot(Improvement,'Color',lines_cmap(i,:),'LineWidth',2)
   lgd2{sub2ind([nNetSizes 1],i,1)} = [num2str(NetSizes(i)),' nodes'];
   hold on
end
xlabel('Iteration/Edges')
ylabel('Average speed up compared to old')
legend(gca,lgd2,'NumColumns',1,'Location','southeast','FontSize',18)
set(gca,'FontSize',18)


print('./images/MatchingDemo3.svg','-dsvg')

clear all

%% Compare the matching generative network model with the old and new code

% Set up a network

Nnodes = 100;
a = zeros(Nnodes);

% Make a pretend distance matrix. We don't need a real one because this is
% a test of speed. 
r = rand(Nnodes);
D = triu(r,1)+triu(r,1)';

% Number of edges to generate. I chose 650 because that is a similar number
% to real world brain networks using a 100 node parcellation in one
% hemisphere
Nedges = 4950;

% Number of iterations to perform
iters = 10;
total_time_new = zeros(iters,1);
total_time_old = zeros(iters,1);

eta = -2;
gam = .4;

iter_time_new = zeros(iters,Nedges);
iter_time_old = zeros(iters,Nedges);

for i = 1:iters

tic
[~,~,iter_time_new(i,:)] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
total_time_new(i) = toc;

tic
[~,~,iter_time_old(i,:)] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
total_time_old(i) = toc;

end

mean_iter_time_new = mean(iter_time_new);
mean_iter_time_old = mean(iter_time_old);

clear i a r D



save('matchingSpeedTestData4.mat')
%

figure('Position',[801 269 1502 559])
subplot(1,3,2)
plot(cumsum(mean_iter_time_new),'LineWidth',2)
hold on
plot(cumsum(mean_iter_time_old),'LineWidth',2)
xlabel('Iteration (network edges)')
ylabel('Cumulative time (seconds)')
legend({'New code','Old code'},'Location','northwest','FontSize',14)
title('Cumulative time')
set(gca,'FontSize',14)
xlim([1 Nedges])

Smoothing_factor = 100;

subplot(1,3,1)
plot(movmean(mean_iter_time_new,Smoothing_factor),'LineWidth',2)
hold on
plot(movmean(mean_iter_time_old,Smoothing_factor),'LineWidth',2)
xlabel('Iteration/Edges')
ylabel('Iteration time (seconds, smoothed)')
legend({'New code','Old code'},'Location','northwest','FontSize',14)
ylimit = ylim;
ylim([ylimit(1) ylimit(2)+((ylimit(2)-ylimit(1))*.15)])
title('Iteration time')
set(gca,'FontSize',14)
xlim([1 Nedges])

subplot(1,3,3)
plot(cumsum(mean_iter_time_old)./cumsum(mean_iter_time_new),'LineWidth',2)
xlabel('Iteration/Edges')
ylabel('Speed up factor of new compared to old code')
title('Speed up factor')
set(gca,'FontSize',14)
xlim([1 Nedges])

print('./images/MatchingDemo4.svg','-dsvg')

clear all

%%

Nnodes = 500;
a = zeros(Nnodes);

% Make a pretend distance matrix. We don't need a real one because this is
% a test of speed. 
r = rand(Nnodes);
D = triu(r,1)+triu(r,1)';

% Number of edges to generate.
Nedges = 124750;

% Number of iterations to perform

eta = -2;
gam = .4;

[~,~,timecourse_new] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);

[~,~,timecourse_old] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);

clear r D a

save('matchingSpeedTestData5.mat')

%

figure('Position',[801 269 1502 559])
subplot(1,3,2)
plot(cumsum(timecourse_new),'LineWidth',2)
hold on
plot(cumsum(timecourse_old),'LineWidth',2)
xlabel('Iteration (network edges)')
ylabel('Cumulative time (seconds)')
legend({'New code','Old code'},'Location','northwest','FontSize',14)
title('Cumulative time')
set(gca,'FontSize',14)
xlim([1 Nedges])

Smoothing_factor = 1000;

subplot(1,3,1)
plot(movmean(timecourse_new,Smoothing_factor),'LineWidth',2)
hold on
plot(movmean(timecourse_old,Smoothing_factor),'LineWidth',2)
xlabel('Iteration/Edges')
ylabel('Iteration time (seconds, smoothed)')
legend({'New code','Old code'},'Location','northwest','FontSize',14)
ylimit = ylim;
ylim([ylimit(1) ylimit(2)+((ylimit(2)-ylimit(1))*.15)])
title('Iteration time')
set(gca,'FontSize',14)
xlim([1 Nedges])

subplot(1,3,3)
plot(cumsum(timecourse_old)./cumsum(timecourse_new),'LineWidth',2)
xlabel('Iteration/Edges')
ylabel('Speed up factor of new compared to old code')
title('Speed up factor')
set(gca,'FontSize',14)
xlim([1 Nedges])

print('./images/MatchingDemo5.svg','-dsvg')

clear all