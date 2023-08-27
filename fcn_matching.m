%function [b,B] = fcn_matching(A,K,D,eta,m,gam,modelvar,epsilon,alpha,AddMult)
load('MatchingTestData.mat')
%AddMult='Mult';

% Function for matching algorithm with comments
% A is an adjacency matrix of seed points (or an empty matrix)
% K is the initial topology matrix
% D is a distance matrix
% m is the number of edges to add
% eta is the parameter for the distance term
% gam is the parameter for the topology term
% modelvar is a 1*2 cell where the first cell contains the model type
% (exponential or powerlaw) for the distance term, the second cell is for
% the topology term
% epsilon is the value to add to the topology term to avoid zeros.
% alpha is the scaling parameter for the topology term. 
% AddMult is either 'Mult' or 'Add' to use that respective model type

% Convert input matrix A to a logical
A = A > 0;

% Add epsilon to matrix K, where K is the inital matching matrix

% If the additive form is used, the epsilon value isn't needed
if strcmp(AddMult,'Add')
    epsilon = 0;
end
K = K + epsilon;

% Get the length of vector Fd. Fd is the distance matrix with the
% powerlaw/exponential decay already applied
n = length(D);

% Calculate the number of edges to add (non-zero elements in A divided by 2)
mseed = nnz(A)/2;

% Get the first element of the input cell array modelvar. This will either
% be 'powerlaw' or 'exponential'
mv1 = modelvar{1};

% Calculate Fd using a power law function or an exponential function using
% the eta parameter
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end

% Get the second element of the input cell array modelvar. This will either
% be 'powerlaw' or 'exponential'
mv2 = modelvar{2};

% Calculate Fk using a power law function or an exponential function using
% the gam parameter
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end

if strcmp(AddMult,'Mult')
    % Calculate Fp (probability function) based on Fd (distance function), Fk
    % (topology function), and complement of A (~A==0 when an edge exists, this
    % just removes that edge from having any impact on the model
    Fp = Fd.*Fk.*~A;
    elseif strcmp(AddMult,'Add')
    % The additive form would look like this I think
        MaxNorm = max(Fk*~A,[],'all');
    if MaxNorm == 0
       TopoTerm = 0; 
    else
       TopoTerm = alpha*(Fk./MaxNorm).*~A;   
    end
       Fp = ((Fd./max(Fd*~A,[],'all')) + TopoTerm).*~A;
end

% Find the upper triangular indices of a ones matrix
[u,v] = find(triu(ones(n),1));

% Calculate indices for an n*n matrix
indx = (v - 1)*n + u;

% Get probabilities P based on Ff and indx (extracts each edge only once)
P = Fp(indx);

% Initialize output vector b
b = zeros(m,1);

% Assign the indices of existing edges in b
b(1:mseed) = find(A(indx));

% Calculate the degree of each node in A
deg = sum(A);

% Create a matrix degmat by replicating deg as columns
degmat = repmat(deg,n,1);

% Transpose degmat
degmat_ = degmat';

% Calculate the sum of degmat and degmat_. This give the summed degree of
% every pair of nodes
degmat_sum = degmat + degmat_;

% The above can also be achieved via the following
% deg = sum(A);
% degmat_sum = (deg+deg').*~eye(n);

% Calculate the neighbor matrix nei by multiplying A with its transpose
nei = (A*A).*~eye(n);

% Loop to find additional edges
for ii = (mseed + 1):m
    
    % Create a cumulative probability distribution
    C = [0; cumsum(P)];
    
    % Select a random value based on the cumulative distribution
    r = sum(rand*C(end) >= C);
    
    % Assign the selected value to b(ii). In this case r represnts the
    % index of the edge to add
    b(ii) = r;
    
    % Get the corresponding node from vector u
    uu = u(r);
    
    % Get the corresponding node from vector v
    vv = v(r);
    
    % Extract the uu-th row of A
    uu_nei = A(uu,:);
    
    % Extract the vv-th row of A
    vv_nei = A(vv,:);
    
    % Connect nodes uu and vv in A
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    
    % Increment the neighbor count for uu and its neighbours
    nei(uu,vv_nei) = nei(uu,vv_nei) + 1;
    nei(vv_nei,uu) = nei(vv_nei,uu) + 1;
    
    % Increment the neighbor count for vv and its neighbours
    nei(vv,uu_nei) = nei(vv,uu_nei) + 1;
    nei(uu_nei,vv) = nei(uu_nei,vv) + 1;
    
    % Increment degmat_sum for uu (rows)
    degmat_sum(uu,:) = degmat_sum(uu,:) + 1;
    
    % Increment degmat_sum for vv (rows)
    degmat_sum(vv,:) = degmat_sum(vv,:) + 1;
    
    % Increment degmat_sum for uu (columns)
    degmat_sum(:,uu) = degmat_sum(:,uu) + 1;
    
    % Increment degmat_sum for vv (columns)
    degmat_sum(:,vv) = degmat_sum(:,vv) + 1;
    
    % This is an easy to follow example of the calculation for the matching
    % index. Basically, double the neighbours for a pair of nodes, then
    % divide this by the sum of the nodes degrees while accounting for any
    % connection that may exist between that node pair
    %K = ( (nei.*2) ./ (degmat_sum - (A.*2) ) ) + epsilon;
    %K(isnan(K)) = epsilon;
    
    % If two nodes have no connections, their matching index will be 0/0
    % which equals nan. We can search and replace nans using 'isnan'
    % however for very large networks, searching for these nans takes a
    % surprising amount of time. To work around this, the section 
    % "(degmat_sum<=2 & nei~=1)" takes the value of 1 when two nodes have 
    % one connection or less each and don't have exactly one neighbor. The
    % reason "degmat_sum<=2" is used is because if two nodes each have a 
    % degree of one but no shared neighbors, this means those two nodes are 
    % connected to each other (and technically share no neighbors). In this 
    % case the equation "degmat_sum - (A.*2)" equals zero (as the summed 
    % degree is cancelled by their shared connection) and could cause the
    % whole equation to fail. The "& nei~=1" catches a case where two nodes
    % are only connected to the same node (i.e., they share one neighbor).
    % If this was not there (i.e., only "degmat_sum<=2" was used) then an  
    % erroneous value of one will be added to the denominator, giving an
    % incorrect result.  
            
%     K = ( (nei.*2) ./ ( (degmat_sum<=2 & nei~=1) + ( degmat_sum - (A.*2) ) ) ) + epsilon;
%     
%     switch mv2
%         case 'powerlaw'
%             Fk = K.^gam;
%         case 'exponential'
%             Fk = exp(gam*K);
%     end

    % Get the nodes to perform the calculation over. Simply find the nodes
    % which the edge was added between (uu vv) and all their neighbours
    % (uu_nei vv_nei). Note I think that uu and vv are not needed here but
    % for safety I have kept them

    % Due to the magic of indexing, we don't actually need to find the
    % unique set of neighbours and can include duplicates
    
    %nodes2use = unique([uu vv find(uu_nei) find(vv_nei)]);
    
    nodes2use = [uu vv find(uu_nei) find(vv_nei)];
    
    % Perform the matching calculation only for the node pairs which can
    % actually be updated.
    
    switch mv2
        case 'powerlaw'
            Fk_update = ( (2 * nei(nodes2use,:) ./ ( (degmat_sum(nodes2use,:)<=2 & nei(nodes2use,:)~=1)+(degmat_sum(nodes2use,:) - (A(nodes2use,:) * 2)) ) ) + epsilon).^gam;
        case 'exponential'
            Fk_update = exp(( (2 * nei(nodes2use,:) ./ ( (degmat_sum(nodes2use,:)<=2 & nei(nodes2use,:)~=1)+(degmat_sum(nodes2use,:) - (A(nodes2use,:) * 2)) ) ) + epsilon)*gam);
    end
        
    % Add in the updated values to the Fk matrix
    Fk(nodes2use,:) = Fk_update;    
    Fk(:,nodes2use) = Fk_update';     
        
    % Recalculate Ff based on updated values
    if strcmp(AddMult,'Mult')
    % Calculate Fp (probability function) based on Fd (distance function), Fk
    % (topology function), and complement of A (~A==0 when an edge exists, this
    % just removes that edge from having any impact on the model
    Fp = Fd.*Fk.*~A;
    elseif strcmp(AddMult,'Add')
    % The additive form would look like this I think
        MaxNorm = max(Fk*~A,[],'all');
    if MaxNorm == 0
       TopoTerm = 0; 
    else
       TopoTerm = alpha*(Fk./MaxNorm);  
    end
        Fp = ((Fd./max(Fd*~A,[],'all')) + TopoTerm).*~A;
    end
        
    % Update probabilities P based on Ff and indx
    P = Fp(indx);
end

% Replace elements of b with their corresponding values from indx
b = indx(b);
B = zeros(n);
B(b) = 1;
B = B + B';