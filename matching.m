function m = matching(A,alt)

% Computes the matching index m for an adjacency matrix A. The matching 
% index is the proportion of a node pairs total connections thats are to
% the same node.

% The alternative way is to define it is out of all the neighbours of a
% pair of nodes, what proportion of them are in common. Set alt = 1 to get
% this definition

if nargin < 2
    alt = 0;
end

% Get the number of nodes n
n = length(A);

% Calculate the neighbours for each pair of nodes
nei = (A*A).*~eye(n);

% Get the degree of each node, followed by the summed degree of each pair
% of nodes
deg = sum(A);

degsum = (deg+deg').*~eye(n);

% Compute the matching index. 
% To avoid dividing by zero we include the term (degsum<=2 & nei~=1). 

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

if alt == 0

m = (nei*2)./( (degsum<=2 & nei~=1) + (degsum-(A.*2))  );

elseif alt == 1

m = (nei)./( (degsum<=2 & nei~=1) + (degsum-(A.*2)-nei)  );
    
end