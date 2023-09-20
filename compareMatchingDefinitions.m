%% COmpare different ways of defining the matching index

% Set the desired number of nodes
n = 1000;

% Make a nested network. See
% https://www.sciencedirect.com/science/article/pii/S037015731930119X for
% deets but in short, a nested network is "a hierarchical organization 
% where the set of neighbors of a node is a subset (superset) of the 
% neighbors of lower (larger) degree". This type of network is also related
% to core-periphery organisation

% Anyways I am using it here because it should sample from the full range
% of possible matching values

A = zeros(n);

for i = 1:n-1
   A(i,1:n+1-i) = 1;
   A(1:n+1-i,i) = 1;
end
A(1:n+1:end) = 0;

% Calculate the neighbours for each pair of nodes
nei = (A*A);

% Get the degree of each node, followed by the summed degree of each pair
% of nodes
deg = sum(A);

degsum = (deg+deg');

% Compute the matching index for the connectivity profiles definition 

m_CP = (nei*2)./( (degsum<=2 & nei~=1) + (degsum-(A.*2))  );

% Compute the matching index for the normalised overlapping neighbourhood profiles definition 

m_NN = (nei)./( (degsum<=2 & nei~=1) + (degsum-(A.*2)-nei)  );

figure

scatter(triu2vec(m_CP,1),triu2vec(m_NN,1),'filled')
xlabel('Matching (connectivity profiles)')
ylabel('Matching (normalised neighbourhood)')
set(gca,'FontSize',14)

print('./images/matchingDefcomparison.svg','-dsvg')

% Convert m_CP to m_NN

CP2NN = (degsum - (A.*2) )./( 2*(degsum - (A.*2) - nei) );

m_NN_cnv = m_CP.*CP2NN;

% Convert m_NN to m_CP

NN2CP = ( 2*(degsum - (A.*2) - nei) )./ (degsum - (A.*2) );

m_CP_cnv = m_NN.*NN2CP;

% Due to precision issues, MATLAB will think m_NN_cnv and m_NN (or m_CP_cnv
% and m_CP) are different results, but I assure you they are exactly the
% same. Break out the pen and paper if you do not believe me!