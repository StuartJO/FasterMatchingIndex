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
        
    r = rand(Nnodes);
    A = triu(r,1)+triu(r,1)';
    A = double(A>(1-den));
    
    % We care about testing speed, so no need to make multiple versions of
    % a network
    
    parfor j = 1:Iters
        tic
        Mold = matching_ind_und(A);
        time_old(i,j) = toc;

        tic
        Mnew = matching(A);
        time_new(i,j) = toc;     
    end
    
    disp(['Completed nodes = ',num2str(Nnodes),', old method = ',num2str(mean(time_old(i,:))),' seconds, new method = ',num2str(mean(time_new(i,:))),' seconds'])
end

figure
subplot(1,3,1)
plot(NetSizes,mean(time_old,2),'LineWidth',4)
hold on
plot(NetSizes,mean(time_new,2),'LineWidth',4)
legend('Old matching method','New matching method')
ylabel('Time in seconds')
xlabel('Number of nodes in network')
set(gca, 'FontSize', 24)
subplot(1,3,2)
plot(NetSizes,mean(time_old,2),'LineWidth',4)
hold on
plot(NetSizes,mean(time_new,2),'LineWidth',4)
legend('Old matching method','New matching method')
ylabel('Time in seconds (log scale)')
xlabel('Number of nodes in network')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 24)
subplot(1,3,3)
plot(NetSizes,mean(time_old,2)./mean(time_new,2),'LineWidth',4)
ylabel('Speed up factor')
xlabel('Number of nodes in network')
set(gca, 'FontSize', 24)

save('matchingSpeedTestData1.mat')
clear all

%% Compare the matching generative network model with the old and new code

% Set up a network

load('random200_data4topomdl.mat')

Nnodes = 100;
seed = zeros(Nnodes);

% Make a pretend distance matrix. We don't need a real one because this is
% a test of speed. 
r = rand(Nnodes);
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

parfor i = 1:iters
tic
B_ = gen_model_mult(seed,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
time_new(i) = toc;
[maxKS_new(i)] = calc_maxKS(A_vals,D,B_); 
end

parfor i = 1:iters
tic
B = gen_model_mult_old(seed,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
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

save('matchingSpeedTestData2.mat')

clear all

%% Compare the matching generative network model with different number of nodes and edges

NetSizes = 100:100:500;
NetEdges = 250:250:1250;

iters = 10;
nNetSizes = length(NetSizes);
nNetEdges = length(NetEdges);

time_new = zeros(nNetSizes,nNetEdges,iters);
time_old = zeros(nNetSizes,nNetEdges,iters);
eta = -2;
gam = .4;

for i = 1:nNetSizes
    
    Nnodes = NetSizes(i);
    a = zeros(Nnodes);
    D = rand(Nnodes);
    D = triu(D,1)+triu(D,1)';

    for j = 1:nNetEdges
        Nedges = NetEdges(j);
        
        parfor k = 1:iters

        tic
        [B,b] = gen_model_mult(a,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
        time_new(i,j,k) = toc;

        tic
        [B_,b_] = gen_model_mult_somewhatold(a,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
        time_old(i,j,k) = toc;
        end  
        
        disp(['Completed nodes = ',num2str(Nnodes),', edges = ',num2str(Nedges),', old method = ',num2str(mean(time_old(i,j,:))),' seconds, new method = ',num2str(mean(time_new(i,j,:))),' seconds'])
        
    end
end

mean_time_new = mean(time_new,3);
mean_time_old = mean(time_old,3);

AverageSpeedUp = mean_time_old./mean_time_new;

figure
imagesc(AverageSpeedUp)
c = colorbar;
xticks(1:5)
yticks(1:5);

xticklabels(NetEdges)
yticklabels(NetSizes)
xlabel('Number of edges')
ylabel('Number of nodes')
c.Label.String = 'Average speed up compared to old';
set(gca,'FontSize',24,'Ydir','normal')

save('matchingSpeedTestData3.mat')

clear all

%% Compare the matching generative network model with the old and new code

NetSizes = [100 250 500 1000 2000];
NetEdges = 500:500:2500;

iters = 10;
nNetSizes = length(NetSizes);
nNetEdges = length(NetEdges);

time_new = zeros(nNetSizes,iters);
time_old = zeros(nNetSizes,iters);

timecourse_new = cell(nNetSizes,1);
timecourse_old = cell(nNetSizes,1);

eta = -2;
gam = .4;

for i = 1:nNetSizes
    
    Nnodes = NetSizes(i);
    a = zeros(Nnodes);
    D = rand(Nnodes);
    D = triu(D,1)+triu(D,1)';
    Nedges = 2500;  
        parfor k = 1:iters

        tStart = tic
        [B,b,t(k,:)] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
        time_new(i,k) = toc(tStart);

        tStart = tic
        [B_,b_,t_(k,:)] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
        time_old(i,k) = toc(tStart);
        end  
        timecourse_new{i} = t;
        timecourse_old{i} = t_;
        disp(['Completed nodes = ',num2str(Nnodes),', edges = ',num2str(Nedges),', old method = ',num2str(mean(time_old(i,:))),' seconds, new method = ',num2str(mean(time_new(i,:))),' seconds'])

end

%

figure
lgd = cell(1,nNetSizes*2);
subplot(1,2,1)

lines_cmap = lines(nNetSizes);

for i = 1:nNetSizes
   newtimes = mean(timecourse_new{i});
   oldtimes = mean(timecourse_old{i});
   plot(cumsum(newtimes),'-','Color',lines_cmap(i,:),'LineWidth',2)
   lgd{sub2ind([2 nNetSizes],1,i)} = ['New code (' num2str(NetSizes(i)),')'];
   hold on
   plot(cumsum(oldtimes),'--','Color',lines_cmap(i,:),'LineWidth',2)
   lgd{sub2ind([2 nNetSizes],2,i)} = ['Old code (' num2str(NetSizes(i)),')'];
end
xlabel('Edges')
ylabel('Cumulative time in seconds')
legend(gca,lgd,'NumColumns',2,'Location','southeast','FontSize',14)
set(gca,'FontSize',18)
set(gca, 'YScale', 'log')

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
xlabel('Edges')
ylabel('Average speed up compared to old')
legend(gca,lgd2,'NumColumns',1,'Location','northeast','FontSize',14)
set(gca,'FontSize',18)

save('matchingSpeedTestData4.mat')

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
time_new = zeros(iters,1);
time_old = zeros(iters,1);

eta = -2;
gam = .4;

t = zeros(iters,Nedges);
t_ = zeros(iters,Nedges);

parfor i = 1:iters

tic
[B,b,t(i,:)] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
time_new(i) = toc;

tic
[B_,b_,t_(i,:)] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
time_old(i) = toc;

end

%disp(['Completed ,',num2str(iters),' networks, old method = ',num2str(mean(time_old)),' seconds, new method = ',num2str(mean(time_new)),' seconds'])
% figure
% yyaxis right
% plot(cumsum(t))
% ylabel('New code')
% 
% yyaxis left
% plot(cumsum(t_))
% ylabel('Old code')

% figure
% 
% plot(cumsum(t,2)','Color','b')
% hold on
% plot(cumsum(t_,2)','Color','r')
% ylabel('Old code')

figure

mean_new_time = mean(t);
mean_old_time = mean(t_);

subplot(1,3,1)
plot(cumsum(mean_new_time))
hold on
plot(cumsum(mean_old_time))

subplot(1,3,2)
plot(mean_new_time)
hold on
plot(mean_old_time)

subplot(1,3,3)
plot(cumsum(mean_old_time)./cumsum(mean_new_time))

save('matchingSpeedTestData5.mat')

clear all

%%

Nnodes = 500;
a = zeros(Nnodes);

% Make a pretend distance matrix. We don't need a real one because this is
% a test of speed. 
r = rand(Nnodes);
D = triu(r,1)+triu(r,1)';

% Number of edges to generate. I chose 650 because that is a similar number
% to real world brain networks using a 100 node parcellation in one
% hemisphere
Nedges = 124750;

% Number of iterations to perform

eta = -2;
gam = .4;

[B,b,timecourse_new] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);

[B_,b_,timecourse_old] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);


figure
subplot(1,3,1)
plot(cumsum(timecourse_new))
hold on
plot(cumsum(timecourse_old))

Smoothing_factor = 1000;

subplot(1,3,2)
plot(movmean(timecourse_new,Smoothing_factor))
hold on
plot(movmean(timecourse_old,Smoothing_factor))

subplot(1,3,3)
plot(cumsum(timecourse_old)./cumsum(timecourse_new))



save('matchingSpeedTestData6.mat')

clear all