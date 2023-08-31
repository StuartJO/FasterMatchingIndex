function [B,b] = gen_model_add(A,PDMs,m,modeltype,modelvar,eta,gam,alpha,normType,epsilon)
% gen_model_add_normmax          Run generative model code for the additive
% model normalising by the max for each term
%
%   Generates synthetic networks using the models described in the study by
%   Oldham et al (2022) in ????.
%
%   Inputs:
%           A,          binary network of seed connections
%           PDMs,       Euclidean distance/fiber length/node similarity
%                       matrix. Multiple can be input either as a cell, 
%                       where each cell contains a different matrix or as a
%                       3D matrix (n*n*nPD, where n is the number of nodes
%                       and nPD is the number of PD matrices).
%           m,          number of connections that should be present in
%                       final synthetic network
%           modeltype,  specifies the generative rule (see below)
%           modelvar,   specifies whether the generative rules are based on
%                       power-law or exponential relationship
%                       ({'powerlaw'}|{'exponential})
%           PDexpo,     the parameter controlling the values in PDMs. If
%                       there are multipe PD matrices, PDexpo should be a
%                       vector where each index gives the marameter for the
%                       corresponding PD matrix
%           gam,        the parameter controlling topology
%           alpha,      the parameter controlling alpha values, should be a
%                       1*3 vector. alpha(1) is the alpha value for the
%                       first PD matrix (should set to 1), alpha(2) is the
%                       value for the topology term, and alpha(3) is for
%                       the second PD matrix
%           normType,   either 'max' or 'sum', this will be the type of
%                       normalisation performed within each respective
%                       term. The default is 'max'
%           epsilon,    the baseline probability of forming a particular
%                       connection (should be a very small number
%                       {default = 0}).
%
%   Output:
%           B,          an adjacency matrix
%           b,          a vector giving the index of each edge in B. Note
%                       that the ordering of b shows which edges formed
%                       first (e.g., b(1) was the fiorst edge to form, b(2)
%                       the second etc etc).
%
%   Full list of model types:
%   (each model type realizes a different generative rule)
%
%       1.  'sptl'          spatial model
%       2.  'neighbors'     number of common neighbors
%       3.  'matching'      matching index
%       4.  'clu-avg'       average clustering coeff.
%       5.  'clu-min'       minimum clustering coeff.
%       6.  'clu-max'       maximum clustering coeff.
%       7.  'clu-diff'      difference in clustering coeff.
%       8.  'clu-prod'      product of clustering coeff.
%       9.  'deg-avg'       average degree
%       10. 'deg-min'       minimum degree
%       11. 'deg-max'       maximum degree
%       12. 'deg-diff'      difference in degree
%       13. 'deg-prod'      product of degree
%       14. 'com'           communicability
%
%       How to convert b to B:
%       n = length(A); B = zeros(n); B(b(:,i)) = 1; B = B + B'; 
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%              Oldham et al (2022) bioRixv
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015
%   Edited by Stuart Oldham, Monash University 2021

if ~exist('epsilon','var')
    epsilon = 0;
end

if ~exist('normType','var')
    normType = 'max';
end
    
if iscell(PDMs) || size(PDMs,3) > 1
    
    if iscell(PDMs)
        if length(PDMs) == 2
             PD = PDMs{2};
             D = PDMs{1};
        else
            PD = [];
            D = PDMs{1};
        end
    else
          PD = squeeze(PDMs(:,:,2));
          D = squeeze(PDMs(:,:,1));
          
    end
else
            PD = [];
            D = PDMs;

end
 n = length(D);
if isempty(PD)
    
    PD = ones(n);
    alpha(3) = 0;
    eta(2) = 1;   
end

alpha(isnan(alpha)) = 0;

switch modeltype

    case 'clu-avg'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@plus,clu(:,ones(1,n)),clu')/2;
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'clu-diff'
        clu = clustering_coef_bu(A);
        Kseed = abs(bsxfun(@minus,clu(:,ones(1,n)),clu'));
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'clu-max'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@max,clu(:,ones(1,n)),clu');
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'clu-min'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@min,clu(:,ones(1,n)),clu');
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'clu-prod'
        clu = clustering_coef_bu(A);
        Kseed = clu*clu';
        b = fcn_clu(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'deg-avg'
        kseed = sum(A,2);
        Kseed = bsxfun(@plus,kseed(:,ones(1,n)),kseed')/2;
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'deg-diff'
        kseed = sum(A,2);
        Kseed = abs(bsxfun(@minus,kseed(:,ones(1,n)),kseed'));
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'deg-max'
        kseed = sum(A,2);
        Kseed = bsxfun(@max,kseed(:,ones(1,n)),kseed');
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'deg-min'
        kseed = sum(A,2);
        Kseed = bsxfun(@min,kseed(:,ones(1,n)),kseed');
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'deg-prod'
        kseed = sum(A,2);
        Kseed = (kseed*kseed').*~eye(n);
        b = fcn_deg(modeltype,A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'neighbors'
        Kseed = (A*A).*~eye(n);
        b = fcn_nghbrs(A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'matching'
        Kseed = matching_ind_und(A);
        b = fcn_matching(A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);

    case 'sptl'
        b = fcn_sptl(A,D,m,eta,modelvar,PD,alpha,normType);

    case 'com'
        Kseed = gexpm(A);
        b = fcn_com(A,Kseed,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType);
end

B = zeros(n);
B(b) = 1;
B = B + B';

function b = fcn_clu(modeltype,A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType)
K = K + epsilon;

n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;

switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end

switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));        
end

switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
switch normType
    case 'max'
        Fk = alpha(2)*(Kf./max(Kf(indx(~EdgesTriu))));
    case 'sum'
        Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
end

Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = (Fd+Fk).*~A;

P = Ff(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    
    switch modeltype
        case 'clu-avg'
    
            K(:,bth) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')/2 + epsilon;
            K(bth,:) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')'/2 + epsilon;

        case 'clu-diff'

            K(:,bth) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)')) + epsilon;
            K(bth,:) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)'))' + epsilon;

        case 'clu-max'
            
            K(:,bth) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
            K(bth,:) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;
            
        case 'clu-min'
            
            K(:,bth) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
            K(bth,:) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;
    
        case 'clu-prod'
            
            K(bth,:) = (c(bth,:)*c') + epsilon;
            K(:,bth) = (c*c(bth,:)') + epsilon;
    end
    
    switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));        
    end
    
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*K);    
    end
    Kf(isinf(Kf)) = 0;
    switch normType
    case 'max'
        Fk = alpha(2)*(Kf./max(Kf(indx(~EdgesTriu))));
    case 'sum'
        Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
    end
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    Ff = (Fd+Fk).*~A;
    P = Ff(indx);
    
end
b = indx(b);


function b = fcn_deg(modeltype,A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType)
n = length(D);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;
Dind = D(indx);
PDind = PD(indx);
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

switch mv1
    case 'powerlaw'
        Df = Dind.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(Dind));    
end
switch mv3
    case 'powerlaw'
        PDf = PDind.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PDind));    
end

switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(~EdgesTriu))) + alpha(3)*(PDf./max(PDf(~EdgesTriu)));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(~EdgesTriu))) + alpha(3)*(PDf./sum(PDf(~EdgesTriu)));        
end


K = K + epsilon;

switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
    switch normType
    case 'max'
        Fk = alpha(2)*(Kf./max(Kf(indx(~EdgesTriu))));
    case 'sum'
        Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
    end
Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    
    
P = (Fd+Fk(indx)).*~A(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    
    
    switch modeltype
        case 'deg-avg'
    
            K(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon];
            K(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon])';
    
        case 'deg-diff'
    
            K(:,w) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon);
            K(w,:) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon)';
    
        case 'deg-min'
            
            K(:,w) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon];
            K(w,:) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]'; 
            
            
        case 'deg-max'
            
            K(:,w) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon];
            K(w,:) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]'; 
            
        case 'deg-prod'
            
            K(:,w) = [k*k(w(1)) + epsilon, k*k(w(2)) + epsilon];
            K(w,:) = [k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]';
    
    end
    
    b(i) = r;
    EdgesTriu(b(i)) = 1;
        
    switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(~EdgesTriu))) + alpha(3)*(PDf./max(PDf(~EdgesTriu)));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(~EdgesTriu))) + alpha(3)*(PDf./sum(PDf(~EdgesTriu)));        
    end
    
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*K);    
    end
    Kf(isinf(Kf)) = 0;
    switch normType
    case 'max'
        Fk = alpha(2)*(Kf./max(Kf(indx(~EdgesTriu))));
    case 'sum'
        Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
    end
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    P = Fd+Fk(indx);

    P(b(1:i)) = 0;
    
end
b = indx(b);

function b = fcn_nghbrs(A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType)
K = K + epsilon;

n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;
switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end

switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));        
end

switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
    switch normType
    case 'max'
        Fk = alpha(2)*(Kf./max(Kf(indx(~EdgesTriu))));
    case 'sum'
        Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
    end
Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
Ff = (Fd+Fk).*~A;

P = Ff(indx);


for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    uu = u(r);
    vv = v(r);
    x = A(uu,:);
    y = A(:,vv);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    K(uu,y) = K(uu,y) + 1;
    K(y,uu) = K(y,uu) + 1;
    K(vv,x) = K(vv,x) + 1;
    K(x,vv) = K(x,vv) + 1;

    switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));        
    end
    
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*K);    
    end
    Kf(isinf(Kf)) = 0;
    
    switch normType
    case 'max'
        Fk = alpha(2)*(Kf./max(Kf(indx(~EdgesTriu))));
    case 'sum'
        Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
    end
    
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    
    Ff = Fd + Fk;

    Ff(A) = 0;

    P = Ff(indx);
    
end
b = indx(b);

function b = fcn_matching(A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType)
% Convert input matrix A to a logical
A = A > 0;

% Add epsilon to matrix K, where K is the inital matching matrix

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
        Df = D.^eta(1);
    case 'exponential'
        Df = exp(eta(1)*D);
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

mv3 = modelvar{2};

switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end

% Find the upper triangular indices of a ones matrix
[u,v] = find(triu(ones(n),1));

% Calculate indices for an n*n matrix
indx = (v - 1)*n + u;

% Initialize output vector b
b = zeros(m,1);

% Assign the indices of existing edges in b
b(1:mseed) = find(A(indx));


EdgesTriu = zeros(size(indx));
EdgesTriu(b(b~=0)) = 1;


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

switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
         MaxNorm = max(Fk*~A,[],'all');
        if MaxNorm == 0
           TopoTerm = 0; 
        else
           TopoTerm = alpha(2)*(Fk./MaxNorm).*~A;   
        end
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu)))); 
        SumNorm = sum(Fk*~A,[],'all')/2;
        if SumNorm == 0
           TopoTerm = 0; 
        else
           TopoTerm = alpha(2)*(Fk./SumNorm).*~A;   
        end
end

    % The additive form would look like this I think
    
       Fp = (Fd + TopoTerm).*~A;

% Get probabilities P based on Ff and indx (extracts each edge only once)
P = Fp(indx);

% Loop to find additional edges
for ii = (mseed + 1):m
    
    % Create a cumulative probability distribution
    C = [0; cumsum(P)];
    
    % Select a random value based on the cumulative distribution
    r = sum(rand*C(end) >= C);
    
    % Assign the selected value to b(ii). In this case r represnts the
    % index of the edge to add
    b(ii) = r;
    
    EdgesTriu(b(ii)) = 1;
    
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
        
    switch normType
        case 'max'
            Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
             MaxNorm = max(Fk*~A,[],'all');
            if MaxNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = alpha(2)*(Fk./MaxNorm).*~A;   
            end
        case 'sum'
            Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu)))); 
            SumNorm = sum(Fk*~A,[],'all')/2;
            if SumNorm == 0
               TopoTerm = 0; 
            else
               TopoTerm = alpha(2)*(Fk./SumNorm).*~A;   
            end
    end
    Fp = (Fd + TopoTerm).*~A;    
    % Update probabilities P based on Ff and indx
    P = Fp(indx);
end
b = indx(b);

function b = fcn_sptl(A,D,m,eta,modelvar,PD,alpha,normType)
n = length(D);
mseed = nnz(A)/2;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

mv1 = modelvar{1};
mv3 = modelvar{3};

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;

switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end
    switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));        
    end

P = Fd(indx).*~A(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;

    EdgesTriu(b(i)) = 1;
    
    switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));        
    end
    
    P = Fd(indx);
    P(b(1:i)) = 0;
end
b = indx(b);

function b = fcn_com(A,K,D,m,eta,gam,modelvar,epsilon,PD,alpha,normType)

K = K + epsilon;

n = length(D);
mseed = nnz(A)/2;
mv1 = modelvar{1};
mv2 = modelvar{2};
mv3 = modelvar{3};

[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;

EdgesTriu = zeros(size(indx));
b = zeros(m,1);
b(1:mseed) = find(A(indx));

EdgesTriu(b(b~=0)) = 1;
switch mv1
    case 'powerlaw'
        Df = D.^eta(1);
    case 'exponential'       
        Df = exp(eta(1)*(D));    
end
switch mv3
    case 'powerlaw'
        PDf = PD.^eta(2);
    case 'exponential'       
        PDf = exp(eta(2)*(PD));    
end
    switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));        
    end
switch mv2
    case 'powerlaw'
        Kf = K.^gam;
    case 'exponential'       
        Kf = exp(gam*K);    
end
Kf(isinf(Kf)) = 0;
    switch normType
    case 'max'
        Fk = alpha(2)*(Kf./max(Kf(indx(~EdgesTriu))));
    case 'sum'
        Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
    end
Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
Ff = (Fd+Fk).*~A;

P = Ff(indx);
%size(P)

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    
    b(i) = r;
    EdgesTriu(b(i)) = 1;
    %r
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    K = gexpm(A)+epsilon;
    switch normType
    case 'max'
        Fd = alpha(1)*(Df./max(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./max(PDf(indx(~EdgesTriu))));
    case 'sum'
        Fd = alpha(1)*(Df./sum(Df(indx(~EdgesTriu)))) + alpha(3)*(PDf./sum(PDf(indx(~EdgesTriu))));        
    end
    switch mv2
        case 'powerlaw'
            Kf = K.^gam;
        case 'exponential'       
            Kf = exp(gam*(K));    
    end
    Kf(isinf(Kf)) = 0;
    switch normType
    case 'max'
        Fk = alpha(2)*(Kf./max(Kf(indx(~EdgesTriu))));
    case 'sum'
        Fk = alpha(2)*(Kf./sum(Kf(indx(~EdgesTriu))));
    end
    Fk(isnan(Fk)) = 0;Fk(isinf(Fk)) = 0;
    Ff = (Fd+Fk).*~A;
    P = Ff(indx);
end
b = indx(b);