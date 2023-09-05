
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

eta = -0.0663620917068697;
gam = 0.268238489046537;

[B_,b_] = gen_model_add(seed,{D},Nedges,'matching',{'exponential','powerlaw'},[-.21,1],[.21;3.96],'max');

[B,b] = gen_model_mult(seed,{D},Nedges,'matching',{'exponential','powerlaw'},eta,gam);

