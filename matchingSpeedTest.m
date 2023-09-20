%% Compare the matching index with different numbers of nodes and a consistent density

% Set the number of iterations
Iters = 10;

% Define a range of network sizes from 100 to 1000 in steps of 100
NetSizes = 100:100:1000;

% Calculate the number of network sizes in the range
nNetSizes = length(NetSizes);

% Create arrays to store execution times for the old and new methods
time_old = zeros(nNetSizes, Iters);
time_new = zeros(nNetSizes, Iters);

% Set the desired density for the random networks
den = 0.2;

% Loop through different network sizes
for i = 1:nNetSizes
    Nnodes = NetSizes(i);

    % Generate a random network with the specified density
    r = rand(Nnodes);
    r_ = triu(r, 1) + triu(r, 1)';
    A = double(r_ > (1 - den));

    % Ensure that the code is loaded into memory by running it once
    if i == 1
        m1 = matching(A);
        m2 = matching_ind_und(A);
    end

    % Perform multiple iterations to measure execution times
    for j = 1:Iters
        tic
        m2 = matching_ind_und(A);
        time_old(i, j) = toc;

        tic
        m1 = matching(A);
        time_new(i, j) = toc;
    end
    
    % Display the average execution times for both methods
    disp(['Completed nodes = ', num2str(Nnodes), ...
        ', old method = ', num2str(mean(time_old(i, :))), ' seconds, new method = ', num2str(mean(time_new(i, :))), ' seconds'])
end

% Clear variables to save memory
clear i j A r r_ Nnodes m1 m2

% Save the data to a file
save('matchingSpeedTestData1.mat')

% Create a figure for plotting
figure('Position', [228 523 1844 563])

% Create subplots for the results
subplot(1, 3, 1)
% Plot the average execution times for both methods
plot(NetSizes, mean(time_new, 2), 'LineWidth', 4)
hold on
plot(NetSizes, mean(time_old, 2), 'LineWidth', 4)
legend('New code', 'Old code', 'Location', 'northwest')
ylabel('Time in seconds')
xlabel('Number of nodes in network')
set(gca, 'FontSize', 18)
xlim([min(NetSizes) max(NetSizes)])

subplot(1, 3, 2)
% Plot the average execution times on a logarithmic scale
plot(NetSizes, mean(time_new, 2), 'LineWidth', 4)
hold on
plot(NetSizes, mean(time_old, 2), 'LineWidth', 4)
legend('New code', 'Old code', 'Location', 'northwest')
ylabel('Time in seconds (log scale)')
xlabel('Number of nodes in network')
set(gca, 'YScale', 'log')
set(gca, 'FontSize', 18)
xlim([min(NetSizes) max(NetSizes)])

subplot(1, 3, 3)
% Plot the speed-up factor of the new code compared to the old code
plot(NetSizes, mean(time_old, 2) ./ mean(time_new, 2), 'LineWidth', 4)
ylabel('Speed up factor of new code')
xlabel('Number of nodes in network')
set(gca, 'FontSize', 18)
xlim([min(NetSizes) max(NetSizes)])

% Save the figure as an SVG file
print('./images/MatchingDemo1.svg', '-dsvg')

% Clear all variables to free up memory
clear all

%% Compare the matching generative network model with the old and new code

% Load data from a file
load('random200_data4topomdl.mat')

% Set the number of nodes in the network
Nnodes = 100;

% Initialize a seed matrix with zeros
seed = zeros(Nnodes);

% Create a pretend distance matrix (for testing purposes)
D = A_dist;

% Define the number of edges to generate (chosen for similarity to real-world brain networks)
Nedges = nnz(adjs{1}) / 2;

% Define the number of iterations to perform
iters = 100;

% Initialize arrays to store execution times for the old and new methods
time_new = zeros(iters, 1);
time_old = zeros(iters, 1);

% Parameters for network generation. These were taken from my 2022 paper
eta = -0.0663620917068697;
gam = 0.268238489046537;

% Extract adjacency matrix and network properties
A = adjs{1};
A_vals{1} = sum(A, 2);
A_vals{2} = clustering_coef_bu(A);
A_vals{3} = betweenness_bin(A)';
A_vals{4} = D(triu(A, 1) > 0);

% Initialize arrays to store maximum KS distances
maxKS_new = zeros(iters, 1);
maxKS_old = zeros(iters, 1);

% Ensure that the code is loaded into memory by running it once
B_ = gen_model_mult(seed, {D}, 1, 'matching', {'exponential', 'powerlaw'}, eta, gam);
B = gen_model_mult_old(seed, {D}, 1, 'matching', {'exponential', 'powerlaw'}, eta, gam);

% Generate networks using both the new and old methods and measure execution times

for i = 1:iters
    tic
    B_ = gen_model_mult(seed, {D}, Nedges, 'matching', {'exponential', 'powerlaw'}, eta, gam);
    time_new(i) = toc;
    [maxKS_new(i)] = calc_maxKS(A_vals, D, B_);
end

for i = 1:iters
    tic
    B = gen_model_mult_old(seed, {D}, Nedges, 'matching', {'exponential', 'powerlaw'}, eta, gam);
    time_old(i) = toc;
    [maxKS_old(i)] = calc_maxKS(A_vals, D, B);
end

% Clear variables to save memory
clear B B_ i adjs

% Save data to a file
save('matchingSpeedTestData2.mat')

% Create a figure for plotting
jitterOffset = 0.5;
figure('Position', [232 262 1040 583])

% Create subplots for time measurements and model fit
subplot(1, 2, 1)

% Scatter plot with jitter for time measurements
cmap = lines(2);
jitter = (rand(length(time_new), 1) - 0.5) * jitterOffset;
x = ones(length(time_new), 1);
scatter(x + jitter, time_new, 'filled', 'MarkerFaceColor', cmap(1, :), 'MarkerFaceAlpha', 0.6);
hold on
jitter = (rand(length(time_old), 1) - 0.5) * jitterOffset;
x = ones(length(time_old), 1) * 2;
scatter(x + jitter, time_old, 'filled', 'MarkerFaceColor', cmap(2, :), 'MarkerFaceAlpha', 0.6);

% Plot mean values with jitter
plotjitterOffset = jitterOffset / 2;
plot([1 - plotjitterOffset 1 + plotjitterOffset], [mean(time_new) mean(time_new)], 'LineWidth', 2, 'Color', 'k')
plot([2 - plotjitterOffset 2 + plotjitterOffset], [mean(time_old) mean(time_old)], 'LineWidth', 2, 'Color', 'k')

% Calculate and display the average speedup
AverageSpeedUp = mean(time_old) / mean(time_new);
title(['Mean speed up = ', num2str(AverageSpeedUp)])
ylabel('Time in seconds')
xticks(1:2)
xticklabels({'New code', 'Old code'})
set(gca, 'FontSize', 18)
xlim([0.5 2.5])

subplot(1, 2, 2)

% Scatter plot with jitter for model fit measurements
cmap = lines(2);
jitter = (rand(length(maxKS_new), 1) - 0.5) * jitterOffset;
x = ones(length(maxKS_new), 1);
scatter(x + jitter, maxKS_new, 'filled', 'MarkerFaceColor', cmap(1, :), 'MarkerFaceAlpha', 0.6);
hold on
jitter = (rand(length(maxKS_old), 1) - 0.5) * jitterOffset;
x = ones(length(maxKS_old), 1) * 2;
scatter(x + jitter, maxKS_old, 'filled', 'MarkerFaceColor', cmap(2, :), 'MarkerFaceAlpha', 0.6);

% Plot mean values with jitter
plotjitterOffset = jitterOffset / 2;
plot([1 - plotjitterOffset 1 + plotjitterOffset], [mean(maxKS_new) mean(maxKS_new)], 'LineWidth', 2, 'Color', 'k')
plot([2 - plotjitterOffset 2 + plotjitterOffset], [mean(maxKS_old) mean(maxKS_old)], 'LineWidth', 2, 'Color', 'k')

xticks(1:2)
xticklabels({'New code', 'Old code'})
ylabel('Model fit')
set(gca, 'FontSize', 18)
xlim([0.5 2.5])

% Save the figure as an SVG file
print('./images/MatchingDemo2.svg', '-dsvg')

% Clear all variables to free up memory
clear all

%% Compare the matching generative network model with the old and new code

% Define different network sizes to test
NetSizes = [100 250 500 1000 2000];

% Number of iterations for each test
iters = 10;
nNetSizes = length(NetSizes);

% Initialize arrays to store execution times
time_new = zeros(nNetSizes, iters);
time_old = zeros(nNetSizes, iters);

% Initialize cell arrays to store execution times for each iteration
timecourse_new = cell(nNetSizes, 1);
timecourse_old = cell(nNetSizes, 1);

% Parameters for network generation
eta = -0.0663620917068697;
gam = 0.268238489046537;

% Create a sample adjacency matrix and distance matrix
a = zeros(10);
D = rand(10);
D = triu(D, 1) + triu(D, 1)';

% Generate networks using both the new and old methods
b1 = matching_gen_model_mult(a, {D}, 1, {'exponential', 'powerlaw'}, eta, gam);
b2 = matching_gen_model_mult_old(a, {D}, 1, {'exponential', 'powerlaw'}, eta, gam);

for i = 1:nNetSizes
    % Set the number of nodes and generate distance matrix
    Nnodes = NetSizes(i);
    a = zeros(Nnodes);
    D = rand(Nnodes);
    D = triu(D, 1) + triu(D, 1)';
    Nedges = 2500;
    
    % Initialize temporary arrays to store execution times for each iteration
    tnew_temp = zeros(iters, Nedges);
    told_temp = zeros(iters, Nedges);
    
    for k = 1:iters
        % Measure execution time for the new method
        tStart = tic;
        [~, ~, tnew_temp(k, :)] = matching_gen_model_mult(a, {D}, Nedges, {'exponential', 'powerlaw'}, eta, gam);
        time_new(i, k) = toc(tStart);
        
        % Measure execution time for the old method
        tStart = tic;
        [~, ~, told_temp(k, :)] = matching_gen_model_mult_old(a, {D}, Nedges, {'exponential', 'powerlaw'}, eta, gam);
        time_old(i, k) = toc(tStart);
    end
    
    % Store execution times for each iteration in cell arrays
    timecourse_new{i} = tnew_temp;
    timecourse_old{i} = told_temp;
    
    % Display information about completed tests
    disp(['Completed nodes = ', num2str(Nnodes), ', edges = ', num2str(Nedges), ...
        ', old method = ', num2str(mean(time_old(i, :))), ' seconds, new method = ', num2str(mean(time_new(i, :))), ' seconds'])
end

% Clear variables to save memory
clear tStart tnew_temp told_temp a D i k b1 b2

% Save data to a file
save('matchingSpeedTestData3.mat')

% Create a figure for plotting
figure('Position', [326 218 1541 652])

subplot(1, 2, 1)

% Create a cumulative time plot for both methods
lines_cmap = lines(nNetSizes);
for i = 1:nNetSizes
    newtimes = mean(timecourse_new{i});
    oldtimes = mean(timecourse_old{i});
    plot(cumsum(newtimes), '-', 'Color', lines_cmap(i, :), 'LineWidth', 2)
    hold on
    plot(cumsum(oldtimes), '--', 'Color', lines_cmap(i, :), 'LineWidth', 2)
end
xlabel('Iteration/Edges')
ylabel('Cumulative time in seconds')
set(gca, 'FontSize', 18)
set(gca, 'YScale', 'log')
xlimits = xlim;
ylimits = ylim;

% Create fake data for plotting a legend
p(1) = plot(-1:-1:-100, -1:-1:-100, '-', 'Color', 'k');
hold on
p(2) = plot(-1:-1:-100, -1:-1:-100, '--', 'Color', 'k');
for i = 1:nNetSizes
    s(i) = scatter(-1, -1, 'filled', 'MarkerEdgeColor', lines_cmap(i, :), 'MarkerFaceColor', lines_cmap(i, :));
    NetSizesLgd{i} = num2str(NetSizes(i));
end
xlim(xlimits)
ylim(ylimits)
lg1 = legend(gca, p, {'New code', 'Old code'}, 'NumColumns', 2, 'Location', 'south', 'FontSize', 18);
ah1 = axes('position', get(gca, 'position'), 'visible', 'off');

lg2 = legend(ah1, s, NetSizesLgd, 'NumColumns', 5, 'Location', 'south', 'FontSize', 18);
xlim(xlimits)
ylim(ylimits)

lg1.Title.String = 'Code version';
lg2.Title.String = 'Network size (nodes)';

lg1_orig_pos = lg1.Position;

lg1.Position(2) = lg1_orig_pos(2) + lg1_orig_pos(4);

subplot(1, 2, 2)
lgd2 = cell(1, nNetSizes);

% Create a plot showing the average speedup
for i = 1:nNetSizes
    newtimes = mean(timecourse_new{i});
    oldtimes = mean(timecourse_old{i});
    Improvement = cumsum(oldtimes) ./ cumsum(newtimes);
    plot(Improvement, 'Color', lines_cmap(i, :), 'LineWidth', 2)
    lgd2{sub2ind([nNetSizes 1], i, 1)} = [num2str(NetSizes(i)), ' nodes'];
    hold on
end
xlabel('Iteration/Edges')
ylabel('Average speed up compared to old')
legend(gca, lgd2, 'NumColumns', 1, 'Location', 'southeast', 'FontSize', 18)
set(gca, 'FontSize', 18)

% Save the figure as an SVG file
print('./images/MatchingDemo3.svg', '-dsvg')

% Clear all variables to free up memory
clear all

%% Compare the model when making all edges of a 100 node network with the old and new code

% Set up a network

% Define the number of nodes in the network
Nnodes = 100;
a = zeros(Nnodes);

% Create a pretend distance matrix. For testing purposes only.
r = rand(Nnodes);
D = triu(r, 1) + triu(r, 1)';

% Define the number of edges to generate. Chosen for similarity to real-world networks.
Nedges = 4950;

% Define the number of iterations to perform
iters = 10;
total_time_new = zeros(iters, 1);
total_time_old = zeros(iters, 1);

% Parameters for network generation
eta = -0.0663620917068697;
gam = 0.268238489046537;

% Initialize arrays to store iteration-specific execution times
iter_time_new = zeros(iters, Nedges);
iter_time_old = zeros(iters, Nedges);

% Generate a network using both the new and old methods (for reference)
b = matching_gen_model_mult(a, {D}, 1, {'exponential', 'powerlaw'}, eta, gam);
b = matching_gen_model_mult_old(a, {D}, Nedges, {'exponential', 'powerlaw'}, eta, gam);

for i = 1:iters
    % Measure execution time for the new method
    tic
    [~, ~, iter_time_new(i, :)] = matching_gen_model_mult(a, {D}, Nedges, {'exponential', 'powerlaw'}, eta, gam);
    total_time_new(i) = toc;
    
    % Measure execution time for the old method
    tic
    [~, ~, iter_time_old(i, :)] = matching_gen_model_mult_old(a, {D}, Nedges, {'exponential', 'powerlaw'}, eta, gam);
    total_time_old(i) = toc;
end

% Calculate the mean iteration time for both methods
mean_iter_time_new = mean(iter_time_new);
mean_iter_time_old = mean(iter_time_old);

% Clear variables to save memory
clear i a r D b

% Save data to a file
save('matchingSpeedTestData4.mat')

% Create a figure for plotting
figure('Position', [801 269 1502 559])

subplot(1, 3, 2)
% Plot cumulative time for both methods
plot(cumsum(mean_iter_time_new), 'LineWidth', 2)
hold on
plot(cumsum(mean_iter_time_old), 'LineWidth', 2)
xlabel('Iteration (network edges)')
ylabel('Cumulative time (seconds)')
legend({'New code', 'Old code'}, 'Location', 'northwest', 'FontSize', 14)
title('Cumulative time')
set(gca, 'FontSize', 14)
xlim([1 Nedges])

% Apply a smoothing factor because individual iterations can be somewhat
% noisy
Smoothing_factor = 100;

subplot(1, 3, 1)
% Plot smoothed iteration time for both methods
plot(movmean(mean_iter_time_new, Smoothing_factor), 'LineWidth', 2)
hold on
plot(movmean(mean_iter_time_old, Smoothing_factor), 'LineWidth', 2)
xlabel('Iteration/Edges')
ylabel('Iteration time (seconds, smoothed)')
legend({'New code', 'Old code'}, 'Location', 'northwest', 'FontSize', 14)
ylimit = ylim;
ylim([ylimit(1) ylimit(2) + ((ylimit(2) - ylimit(1)) * 0.15)])
title('Iteration time')
set(gca, 'FontSize', 14)
xlim([1 Nedges])

subplot(1, 3, 3)
% Plot the speedup factor of the new code compared to the old code
plot(cumsum(mean_iter_time_old) ./ cumsum(mean_iter_time_new), 'LineWidth', 2)
xlabel('Iteration/Edges')
ylabel('Speed up factor of new compared to old code')
title('Speed up factor')
set(gca, 'FontSize', 14)
xlim([1 Nedges])

% Save the figure as an SVG file
print('./images/MatchingDemo4.svg', '-dsvg')

% Clear all variables to free up memory
clear all

%% Compare the model when making all edges of a 500 node network with the old and new code

% Set the number of nodes in the network
Nnodes = 500;
a = zeros(Nnodes);

% Create a pretend distance matrix. For testing purposes only.
r = rand(Nnodes);
D = triu(r, 1) + triu(r, 1)';

% Define the number of edges to generate.
Nedges = 124750;

% Parameters for network generation
eta = -0.0663620917068697;
gam = 0.268238489046537;

% Generate a network using both the new and old methods (for reference)
b = matching_gen_model_mult(a, {D}, 1, {'exponential', 'powerlaw'}, eta, gam);
b = matching_gen_model_mult_old(a, {D}, Nedges, {'exponential', 'powerlaw'}, eta, gam);

% Measure execution time for the new method
[~, ~, timecourse_new] = matching_gen_model_mult(a, {D}, Nedges, {'exponential', 'powerlaw'}, eta, gam);

% Measure execution time for the old method
[~, ~, timecourse_old] = matching_gen_model_mult_old(a, {D}, Nedges, {'exponential', 'powerlaw'}, eta, gam);

% Clear variables to save memory
clear r D a b 

% Save data to a file
save('matchingSpeedTestData5.mat')

% Create a figure for plotting
figure('Position', [801 269 1502 559])
subplot(1, 3, 2)
% Plot cumulative time for both methods
plot(cumsum(timecourse_new), 'LineWidth', 2)
hold on
plot(cumsum(timecourse_old), 'LineWidth', 2)
xlabel('Iteration (network edges)')
ylabel('Cumulative time (seconds)')
legend({'New code', 'Old code'}, 'Location', 'northwest', 'FontSize', 14)
title('Cumulative time')
set(gca, 'FontSize', 14)
xlim([1 Nedges])

Smoothing_factor = 1000;

subplot(1, 3, 1)
% Plot smoothed iteration time for both methods
plot(movmean(timecourse_new, Smoothing_factor), 'LineWidth', 2)
hold on
plot(movmean(timecourse_old, Smoothing_factor), 'LineWidth', 2)
xlabel('Iteration/Edges')
ylabel('Iteration time (seconds, smoothed)')
legend({'New code', 'Old code'}, 'Location', 'northwest', 'FontSize', 14)
ylimit = ylim;
ylim([ylimit(1) ylimit(2) + ((ylimit(2) - ylimit(1)) * 0.15)])
title('Iteration time')
set(gca, 'FontSize', 14)
xlim([1 Nedges])

subplot(1, 3, 3)
% Plot the speedup factor of the new code compared to the old code
plot(cumsum(timecourse_old) ./ cumsum(timecourse_new), 'LineWidth', 2)
xlabel('Iteration/Edges')
ylabel('Speed up factor of new compared to old code')
title('Speed up factor')
set(gca, 'FontSize', 14)
xlim([1 Nedges])

% Save the figure as an SVG file
print('./images/MatchingDemo5.svg', '-dsvg')

% Clear all variables to free up memory
clear all
