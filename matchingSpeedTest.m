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

Nnodes = 100;
a = zeros(Nnodes);

% Make a pretend distance matrix. We don't need a real one because this is
% a test of speed. 
r = rand(Nnodes);
D = triu(r,1)+triu(r,1)';

% Number of edges to generate. I chose 650 because that is a similar number
% to real world brain networks using a 100 node parcellation in one
% hemisphere
Nedges = 650;

% Number of iterations to perform
iters = 100;
time_new = zeros(iters,1);
time_old = zeros(iters,1);

eta = -2;
gam = .4;

parfor i = 1:iters

tic
[B,b] = gen_model_mult(a,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
time_new(i) = toc;

tic
[B_,b_] = gen_model_mult_old(a,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
time_old(i) = toc;

end

disp(['Completed ,',num2str(iters),' networks, old method = ',num2str(mean(time_old)),' seconds, new method = ',num2str(mean(time_new)),' seconds'])
        

figure
boxplot([time_new time_old])
AverageSpeedUp = mean(time_old)/mean(time_new);
title(['Mean speed up = ',num2str(AverageSpeedUp)])
ylabel('Time in seconds')

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
        [B_,b_] = gen_model_mult_old(a,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);
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
iters = 1;
time_new = zeros(iters,1);
time_old = zeros(iters,1);

eta = -2;
gam = .4;

for i = 1:iters

tic
[B,b,t] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
time_new(i) = toc;

tic
[B_,b_,t_] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
time_old(i) = toc;

end

%disp(['Completed ,',num2str(iters),' networks, old method = ',num2str(mean(time_old)),' seconds, new method = ',num2str(mean(time_new)),' seconds'])
figure
yyaxis right
plot(cumsum(t))
ylabel('New code')

yyaxis left
plot(cumsum(t_))
ylabel('Old code')

figure
boxplot([time_new time_old])
AverageSpeedUp = mean(time_old)/mean(time_new);
title(['Mean speed up = ',num2str(AverageSpeedUp)])
ylabel('Time in seconds')

save('matchingSpeedTestData5.mat')

clear all

%% Compare the matching generative network model with the old and new code

NetSizes = 100:100:500;
NetEdges = 250:250:1250;

iters = 10;
nNetSizes = length(NetSizes);
nNetEdges = length(NetEdges);

time_new = zeros(nNetSizes,nNetEdges);
time_old = zeros(nNetSizes,nNetEdges);

timecourse_new = cell(nNetSizes,nNetEdges);
timecourse_old = cell(nNetSizes,nNetEdges);

eta = -2;
gam = .4;

for i = 1:nNetSizes
    
    Nnodes = NetSizes(i);
    a = zeros(Nnodes);
    D = rand(Nnodes);
    D = triu(D,1)+triu(D,1)';

    for j = 1:nNetEdges
        Nedges = NetEdges(j);
        
        t = zeros(iters,Nedges);
        t_ = zeros(iters,Nedges);
      
        parfor k = 1:iters

        tic
        [B,b,t(k,:)] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
        time_new(i,j,k) = toc;

        tic
        [B_,b_,t_(k,:)] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
        time_old(i,j,k) = toc;
        end  
        timecourse_new{i,j} = t;
        timecourse_old{i,j} = t_;
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

%

figure
for i = 1:nNetSizes
    for j = 1:nNetEdges
        subplot(nNetSizes,nNetEdges,sub2ind([nNetSizes nNetEdges],j,i))
        newtimes = mean(timecourse_new{i,j});
        oldtimes = mean(timecourse_old{i,j});
        plot(cumsum(newtimes)); hold on; plot(cumsum(oldtimes))
        xlabel('Iteration')
        ylabel('Cumulative time in seconds')
        legend('New code','Old code','Location','northwest')
        title([num2str(NetSizes(i)),' nodes, ',num2str(NetEdges(j)),' edges'])
        xlim([1 NetEdges(j)])
        
        mean_total_time_new(i,j) = sum(newtimes);
        mean_total_time_old(i,j) = sum(oldtimes);
        
        ylim([0 max(mean_total_time_new(i,j),mean_total_time_old(i,j))])        
    end
end

AverageSpeedUpTotal = mean_total_time_old./mean_total_time_new;

figure
imagesc(AverageSpeedUpTotal)
c = colorbar;
xticks(1:5)
yticks(1:5);

xticklabels(NetEdges)
yticklabels(NetSizes)
xlabel('Number of edges')
ylabel('Number of nodes')
c.Label.String = 'Average speed up compared to old';
set(gca,'FontSize',24,'Ydir','normal')

save('matchingSpeedTestData6.mat')

clear all



%%
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

        tic
        [B,b,t(k,:)] = matching_gen_model_mult(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
        time_new(i,k) = toc;

        tic
        [B_,b_,t_(k,:)] = matching_gen_model_mult_old(a,{D},Nedges,{'exponential','powerlaw'},eta,gam);
        time_old(i,k) = toc;
        end  
        timecourse_new{i} = t;
        timecourse_old{i} = t_;
        disp(['Completed nodes = ',num2str(Nnodes),', edges = ',num2str(Nedges),', old method = ',num2str(mean(time_old(i,:))),' seconds, new method = ',num2str(mean(time_new(i,:))),' seconds'])

end

%

figure
lgd = cell(1,nNetSizes*2);
subplot(1,2,1)
for i = 1:nNetSizes
    newtimes = mean(timecourse_new{i});
   oldtimes = mean(timecourse_old{i});
   plot(cumsum(newtimes))
   lgd{sub2ind([nNetSizes 2],i,1)} = ['New code,' num2str(NetSizes(i)),' nodes'];
   hold on
   plot(cumsum(newtimes),'--')
   hold on
   lgd{sub2ind([nNetSizes 2],i,2)} = ['Old code,' num2str(NetSizes(i)),' nodes'];
end
xlabel('Edges')
ylabel('Cumulative time in seconds')
legend(ldg,'NumColumns',nNetSizes,'Location','northwest')

subplot(1,2,1)
lgd = cell(1,nNetSizes);
subplot(1,2,1)
for i = 1:nNetSizes
    newtimes = mean(timecourse_new{i});
    oldtimes = mean(timecourse_old{i});
   Improvement = oldtimes./newtimes;
   plot(Improvement)
   lgd{sub2ind([nNetSizes 1],i,1)} = [num2str(NetSizes(i)),' nodes'];
   hold on
end
xlabel('Edges')
ylabel('Average speed up compared to old')
legend(ldg,'Location','northwest')


save('matchingSpeedTestData6.mat')

clear all
